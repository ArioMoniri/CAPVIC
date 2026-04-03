"""Shared async HTTP client infrastructure with rate limiting and retry logic."""

from __future__ import annotations

import asyncio
import logging
import sys
from typing import Any

import httpx

from variant_mcp.constants import MAX_RETRIES, REQUEST_TIMEOUT, RETRY_BACKOFF_BASE

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stderr))


class RateLimiter:
    """Token-bucket style rate limiter using asyncio.Semaphore."""

    def __init__(self, max_per_second: int) -> None:
        self._semaphore = asyncio.Semaphore(max_per_second)
        self._delay = 1.0 / max_per_second

    async def acquire(self) -> None:
        await self._semaphore.acquire()
        asyncio.get_running_loop().call_later(self._delay, self._semaphore.release)


class BaseClient:
    """Base async HTTP client with rate limiting, retries, and error handling."""

    def __init__(
        self,
        base_url: str,
        rate_limit: int,
        timeout: int = REQUEST_TIMEOUT,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self._rate_limiter = RateLimiter(rate_limit)
        self._timeout = timeout
        self._headers = headers or {}
        self._client: httpx.AsyncClient | None = None

    async def _get_client(self) -> httpx.AsyncClient:
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(self._timeout),
                headers=self._headers,
                follow_redirects=True,
            )
        return self._client

    async def close(self) -> None:
        if self._client and not self._client.is_closed:
            await self._client.aclose()

    async def _request(
        self,
        method: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        json_body: dict[str, Any] | None = None,
        max_retries: int = MAX_RETRIES,
    ) -> httpx.Response:
        """Execute HTTP request with rate limiting and exponential backoff retry."""
        client = await self._get_client()
        last_error: Exception | None = None

        for attempt in range(max_retries + 1):
            await self._rate_limiter.acquire()
            try:
                response = await client.request(
                    method,
                    url,
                    params=params,
                    json=json_body,
                )
                response.raise_for_status()
                return response
            except httpx.HTTPStatusError as e:
                last_error = e
                if e.response.status_code in (429, 500, 502, 503, 504):
                    wait = RETRY_BACKOFF_BASE * (2**attempt)
                    logger.warning(
                        "HTTP %d from %s, retrying in %.1fs (attempt %d/%d)",
                        e.response.status_code,
                        url,
                        wait,
                        attempt + 1,
                        max_retries,
                    )
                    await asyncio.sleep(wait)
                    continue
                raise ClientError(
                    f"HTTP {e.response.status_code} from {url}: {e.response.text[:500]}"
                ) from e
            except httpx.TimeoutException as e:
                last_error = e
                if attempt < max_retries:
                    wait = RETRY_BACKOFF_BASE * (2**attempt)
                    logger.warning("Timeout for %s, retrying in %.1fs", url, wait)
                    await asyncio.sleep(wait)
                    continue
            except httpx.RequestError as e:
                last_error = e
                if attempt < max_retries:
                    wait = RETRY_BACKOFF_BASE * (2**attempt)
                    logger.warning("Request error for %s: %s, retrying", url, e)
                    await asyncio.sleep(wait)
                    continue

        raise ClientError(
            f"Failed after {max_retries + 1} attempts to {url}: {last_error}"
        ) from last_error

    async def get(self, path: str, params: dict[str, Any] | None = None) -> httpx.Response:
        url = f"{self.base_url}/{path.lstrip('/')}" if path else self.base_url
        return await self._request("GET", url, params=params)

    async def post(
        self,
        path: str,
        json_body: dict[str, Any] | None = None,
        params: dict[str, Any] | None = None,
    ) -> httpx.Response:
        url = f"{self.base_url}/{path.lstrip('/')}" if path else self.base_url
        return await self._request("POST", url, params=params, json_body=json_body)

    async def get_json(self, path: str, params: dict[str, Any] | None = None) -> Any:
        response = await self.get(path, params=params)
        return response.json()

    async def post_json(self, path: str, json_body: dict[str, Any] | None = None) -> Any:
        url = f"{self.base_url}/{path.lstrip('/')}" if path else self.base_url
        response = await self._request("POST", url, json_body=json_body)
        return response.json()


class ClientError(Exception):
    """Raised when an API client encounters an error."""
