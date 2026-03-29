# F3 Falsification Report — Background Response

- Hard fail: **YES**
- Included families: ['flrw', 'schwarzschild']
- Rule: fail_ratio >= 0.5, top-k requirement=2
- N grid: [256, 512]

## Hard-fail items
- flrw param=1.0: fail_ratio=0.50 (1/2)

## Notes

- This is a skeleton runner; metric-faithful background generators can be swapped in-place.
- Current script explicitly keeps weak/moderate windows pre-registered in config.