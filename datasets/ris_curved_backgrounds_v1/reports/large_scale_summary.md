# Large-Scale Curvature Background Robustness Report

**Generated:** 2026-03-28 13:39:25  
**Scope:** N=512-2048 (3 partitions actually executed: 512/1024/2048)  
**Total Trials:** 180 (3 sizes × 6 configurations × 10 trials)

## Summary Statistics

### Partition N=512
- All 6 configurations success rate: **100.0%**
- Rank diversity: 10 distinct nodes per configuration

### Partition N=1024
- All 6 configurations success rate: **100.0%**
- Rank diversity: 10 distinct nodes per configuration

### Partition N=2048
- All 6 configurations success rate: **100.0%**
- Rank diversity: 10 distinct nodes per configuration

## Key Findings

1. Mild curvature (FLRW κ≤1.0, Schwarzschild φ≤0.01) ranking stable.
2. Moderate/strong settings in this run still preserved 100% top-node-in-seeds.
3. Scaling to N=2048 did not degrade success rate.
