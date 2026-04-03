# New Implementation Checklist

When adding a new data source, API client, or MCP tool to CAPVIC, follow this checklist to ensure all files are updated and all tests pass.

---

## 1. API Client

- [ ] Create `src/variant_mcp/clients/{source}_client.py`
  - Extend `BaseClient` with `base_url`, `rate_limit`, `headers`
  - Implement async methods for each API endpoint
  - Add proper error handling (try/except with `logger.warning`, return empty/None on failure)
  - Handle API quirks (POST vs GET, URL encoding, trailing slashes, response type coercion)
- [ ] Add constants to `src/variant_mcp/constants.py`
  - `{SOURCE}_API_URL` — base URL
  - `{SOURCE}_RATE_LIMIT` — requests per second (be courteous)
  - Any additional URLs (e.g., v1 vs v2 endpoints)
- [ ] Register in `src/variant_mcp/clients/__init__.py`
  - Add import
  - Add to `__all__` list

## 2. Data Models

- [ ] Add Pydantic models to `src/variant_mcp/models/evidence.py`
  - Model for the API response data (e.g., `CancerHotspot`, `DriverMutationAssessment`)
  - All fields should have sensible defaults (Optional, default=None, default_factory=list)
- [ ] Update `EvidenceBundle` model if the data feeds into classification
  - Add field (e.g., `cancer_hotspots: list[CancerHotspot] = Field(default_factory=list)`)
  - Add computed property (e.g., `has_hotspot_data`)
  - Add to `sources_queried` property

## 3. MCP Tool Registration

- [ ] Add tool function(s) to `src/variant_mcp/server.py`
  - Use `@mcp.tool()` decorator with rich description including "Use for:" examples
  - Accept `output_format` parameter (markdown/json/text)
  - Format response based on output_format
  - Handle edge cases (empty results, API failures)
- [ ] Create singleton client instance at module level (e.g., `litvar_client = LitVarClient()`)
- [ ] Update server instructions string (the system prompt that lists available tools)

## 4. Evidence Pipeline Integration (if applicable)

- [ ] Wire into `_gather_evidence()` in `server.py`
  - Add async task to the `tasks` dict
  - Add result parsing in the evidence bundle loop
- [ ] Wire into classification scorers if relevant
  - Update `classification/oncogenicity_sop.py` if the data informs oncogenicity codes
  - Update `classification/amp_asco_cap.py` if it informs AMP tiers

## 5. Tests

- [ ] Create or extend test file in `tests/`
  - Mock all HTTP calls with `respx` — never hit real APIs in tests
  - Test the client methods (happy path, empty result, API error)
  - Test data model creation and field access
  - Test evidence bundle integration (if applicable)
  - Test classification impact (if applicable)
  - **Watch for URL encoding issues** — verify mock URLs match what httpx actually sends
  - **Watch for trailing slashes** — some APIs return 404 with/without them
- [ ] Run full test suite: `.venv/bin/python -m pytest tests/ --tb=short`
- [ ] Run lint: `.venv/bin/ruff check src/ tests/`
- [ ] Run format check: `.venv/bin/ruff format --check src/ tests/`

## 6. Documentation Updates

### README.md

- [ ] Update **Overview > What It Does** table (add capability row)
- [ ] Update **Data Sources** table (add source row with access info)
- [ ] Update **Data Sources** SVG alt text (e.g., "8 Integrated" -> "10 Integrated")
- [ ] Update **Project Structure** tree (add client file)
- [ ] Update tool count references (search for old count, e.g., "26 tools")
- [ ] Update **Tool Reference** section (add new category or add to existing)
- [ ] Update **Natural Language Prompts > Example Prompts** table
- [ ] Update **Example Workflows** section (add new workflow with real output)
- [ ] Update **Data Sources & Freshness** table
- [ ] Update **Genome Build Reference** table (if build-sensitive)
- [ ] Update **Bioinformatician's Assessment** numbered list
- [ ] Update **Key References** table (add paper citations)
- [ ] Update **OpenCode** section tool count
- [ ] Update **Docker > Example Test Commands** table (add JSON-RPC example)

### CHANGELOG.md

- [ ] Add version section with Added/Changed/Fixed subsections
- [ ] List all new clients, models, tools, tests
- [ ] Note breaking changes (e.g., type widening in base_client)
- [ ] Note bug fixes discovered during implementation

### ROADMAP.md

- [ ] Move items from "Planned" to "Completed" under the new version
- [ ] Add future enhancement ideas to the appropriate planned version

### Other docs

- [ ] Update `docs/example-prompts-real-answers.md` with real API responses
- [ ] Verify all example outputs are current (APIs change over time)

## 7. SVG Assets

All SVGs in `assets/` must be regenerated when sources or tools change:

- [ ] `assets/architecture.svg` — update tool count, add module boxes, add data source cards, update version
- [ ] `assets/data-sources.svg` — add source cards, update count in title and footer
- [ ] `assets/data-flow.svg` — add evidence bundle rows, update pipeline steps, update version
- [ ] `assets/tools-overview.svg` — add tool entries to categories, update counts, update version
- [ ] `assets/classification-frameworks.svg` — update only if new framework added
- [ ] Update version string in ALL SVG footers (e.g., "v1.0.5" -> "v1.0.6")

## 8. Version Bump

- [ ] Update version in `pyproject.toml` (`version = "X.Y.Z"`)
- [ ] Verify version matches in all SVG footers

## 9. Pre-Commit Verification

```bash
# All must pass before committing:
.venv/bin/ruff check src/ tests/
.venv/bin/ruff format --check src/ tests/
.venv/bin/python -m pytest tests/ --tb=short
python -c "from variant_mcp.server import mcp; print('Server OK')"
```

## 10. Live API Verification

Before releasing, test against real APIs (not just mocks):

```bash
.venv/bin/python -c "
import asyncio
from variant_mcp.clients.{source}_client import {Source}Client

async def test():
    c = {Source}Client()
    result = await c.{method}('{test_input}')
    print(result)

asyncio.run(test())
"
```

- [ ] Verify API returns expected structure
- [ ] Check for URL encoding issues (especially `#`, `%`, `@` in IDs)
- [ ] Check trailing slash sensitivity
- [ ] Check POST vs GET method
- [ ] Check JSON body format (array vs object)
- [ ] Check error handling for nonexistent inputs
- [ ] Check rate limiting behavior
- [ ] Document any API quirks in the example prompts file

## 11. CI/CD

- [ ] Commit with descriptive message
- [ ] Push and verify CI passes on both Python 3.11 and 3.12
- [ ] If CI fails, check the test that uses the new mock — URL encoding mismatches are the #1 cause

---

## Common Pitfalls

| Pitfall | Example | Fix |
|---------|---------|-----|
| `#` in URLs | LitVar IDs contain `##` — httpx strips as fragment | Pre-encode as `%23` before building URL path |
| Trailing slash | LitVar2 `/variant/get/{id}/` returns 404 | Remove trailing slash from path |
| POST not GET | Cancer Hotspots requires POST with JSON array body | Use `self.post()` with `json_body=[gene]` |
| json_body type | Cancer Hotspots sends `["GENE"]` not `{"gene":"GENE"}` | Widen `json_body` type to `Any` in base_client |
| q-value strings | Cancer Hotspots returns `"3.65e-82"` as string | `float(raw_q) if raw_q is not None else None` |
| Mock URL mismatch | respx mock URL must match what httpx actually sends | Print the actual URL httpx sends and use that in mock |
| Empty/null fields | API returns `null` or `{}` for optional fields | Use `or {}` / `or []` defensive patterns |
| API version drift | Endpoints change (e.g., ClinVar deprecated `rettype=variation`) | Test against live API periodically |
