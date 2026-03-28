# Cost-Benefit Analysis: N>512 Curvature Background Testing

## Parameter Breakdown

- N=512: 2.6 MB, ~2.6s (10 trials)
- N=768: 5.8 MB, ~5.9s
- N=1024: 10.4 MB, ~10.5s
- N=1536: 23.4 MB, ~23.6s
- N=2048: 41.6 MB, ~41.9s

## Recommendations

- Recommended max in current setup: **N=2048**
- Suggested batch size: 10-15 configs/partition
- Keep >30% free RAM buffer
