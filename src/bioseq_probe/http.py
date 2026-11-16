"""A polite HTTP client for the public EBI/UniProt APIs, with retries."""

import time

import httpx

CLIENT = httpx.Client(timeout=120, headers={"User-Agent": "llm-bioseq-probe/0.1"})


def get(url: str, params: dict | None = None, tries: int = 10) -> httpx.Response:
    for attempt in range(tries):
        time.sleep(0.2)  # stay under the EBI rate limits
        try:
            r = CLIENT.get(url, params=params)
            if r.status_code == 429:
                time.sleep(float(r.headers.get("retry-after", 30)))
            if r.status_code in (408, 429) or r.status_code >= 500:
                raise httpx.HTTPStatusError(str(r.status_code), request=r.request, response=r)
            r.raise_for_status()
            return r
        except (httpx.TransportError, httpx.HTTPStatusError):
            if attempt == tries - 1:
                raise
            time.sleep(min(2**attempt, 60))
    raise AssertionError
