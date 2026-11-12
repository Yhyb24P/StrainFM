# Expanded AMR joint split (pre-registered)

Joint near-clone-isolated split over identity representatives labeled from the BV-BRC AST snapshot. Whole near-clone groups form components; source components (declared projects, ENA studies, BioProject) are reported per partition but not isolated, because the declared-source graph is one connected blob. The same BioSample keeps one partition across drugs. Frozen before any model result is inspected.

- ciprofloxacin: N=3218 (S=1884, R=1334), groups=1970, components=40, conflicts=24 groups; partitions train: 1931 primary units / 1376 groups / 39 source components (largest 0.453) / S1291 R956; validation: 286 primary units / 200 groups / 24 source components (largest 0.47) / S188 R140; test: 620 primary units / 394 groups / 29 source components (largest 0.467) / S405 R238
- gentamicin: N=3272 (S=2464, R=808), groups=2054, components=39, conflicts=51 groups; partitions train: 1912 primary units / 1439 groups / 38 source components (largest 0.446) / S1689 R607; validation: 272 primary units / 198 groups / 22 source components (largest 0.483) / S253 R68; test: 609 primary units / 417 groups / 27 source components (largest 0.46) / S522 R133

Old AMR-02 artifacts are unchanged; overlap with the seen AMR-02 test must be reported in any confirmatory claim.

## 验收检查（ds，2026-11-12）

- 跨药一致性：3,456 个基因组在 cipro/gent 两任务中的分区完全一致（0 不一致）。
- 来源分层：各分区来源组件 22–39 个，最大单来源占比 44.6%–48.3%；来源组件跨分区，**不支持“来源独立测试”主张**。
- 与旧 AMR-02 已见测试重叠：旧 test 193 个 BioSample 中 22 个出现在新队列，其中 **3 个进入新 test**（约 0.5%）；确认性结论须据此限定，不静默删改已冻结分区。
- 严格来源隔离的退化结果（23 组件、val 24 例）保留在 `../escherichia_coli_ast_joint_split_source_isolated_20261112/`。
