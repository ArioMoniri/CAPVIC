FROM python:3.11-slim

WORKDIR /app

# Copy source first, then install (hatchling needs src/ to build)
COPY pyproject.toml .
COPY src/ src/

RUN pip install --no-cache-dir .

EXPOSE 8000

# MCP stdio transport — the server reads JSON-RPC from stdin, writes to stdout
CMD ["variant-mcp"]
