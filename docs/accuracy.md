# 精度约定与容差

下列容差在 `orbit_service/constants.py` 中定义，并由自动化测试守护。

| 项目 | 约定容差 | 测试 |
| --- | --- | --- |
| ISO 字符串与 Unix 时间戳一致性 | < 1 µs | `test_inputs` |
| TLE 根数与 COE 描述同一轨道 | < 1 m（换算层面） | `test_tle_and_coe_equivalent` |
| J2 外推 vs 独立 RK4 参考轨道（近地，2 天，默认 15 s 步长） | < 50 m | `test_j2_matches_independent_reference` |
| 无摄动模式 7 天比能量相对漂移 | < 1×10⁻⁸ | `test_two_body_energy_and_momentum_conserved` |
| 无摄动模式 7 天角动量模相对漂移 | < 1×10⁻¹⁰ | 同上 |
| 窗口边界收紧后稳定性 | 相邻容差级内移动有界，窗口不跳变/丢失/乱序 | `test_windows.py` |

说明：

- “独立参考轨道”指测试中用**独立书写的 RK4 循环**、以 1 s 步长积分同一
  J2 动力学生成的轨迹（见 `tests/conftest.py`）。它与生产代码共用加速度
  定义但不共用积分实现，用于检验积分收敛与步长选择。
- 50 m 容差针对默认 15 s 步长；缩短 `step_seconds` 可按四阶规律显著提高
  精度。
- TLE 与 COE 的“一致”是根数换算层面的一致（n → a，角度换算），两者汇入
  同一内部平均根数；不含真实 SGP4 与不同大气/光压模型之间的差异。
