# Data governance

Real TRACERx data on `/Nemo` is patient data under strict access controls. Non-negotiable rules.

## Development

- Local dev uses `nemo_mock/` only. Synthetic, non-identifiable.
- The agent never sees real data during code development.
- Real-data execution happens exclusively on the HPC with `TX_DATA_ROOT=/Nemo`.

## LLM boundary (when we add the assistant)

- The LLM never receives raw patient rows.
- The LLM sees only: catalog YAML (schema + semantics), summary statistics, and aggregated query results.
- LLM-emitted code executes in a sandbox against real data; only aggregated / structured results return.
- API-based LLMs (Claude, GPT-4, Gemini) may only be used after data-access-committee review.
- Local LLMs (Llama-3.1 on HPC) preferred for any workflow touching identifiable data.

## Repository hygiene

- Never commit real data, real patient IDs, or real filenames from the cohort.
- `nemo_mock/` uses synthetic IDs (`patient_mock_001`, etc.) and synthetic values.
- If a real path or ID slips into a commit, stop and flag to User before pushing anywhere.

## Provenance

Every derived artefact records what it was built from (source SHA256) and by what code (repo-relative script path + commit). Traceability is required for anything that touches real data.
