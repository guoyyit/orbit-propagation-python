# Local Orbit Propagation & Ground-Station Visibility Service

纯本地运行的卫星轨道预报与地面站可见性计算服务（FastAPI）。支持 TLE 与
经典轨道六根数输入、J2 摄动外推、可见窗口、地影事件、批量并发、确定性
复现与场景持久化。**不依赖任何外部服务，结果可复现。**

## 快速开始

```bash
uv sync
uv run uvicorn app.main:app --reload
# 测试
uv run pytest -q
```

需要 Python 3.14（见 `pyproject.toml`）。

## 功能

* **输入**：TLE（严格格式/校验和/字段范围校验，可区分错误码与位置）或
  经典六根数（角度弧度，SI 单位）。
* **外推**：任意时刻位置/速度；时间支持 ISO 8601 与 Unix 秒（结果一致）；
  二体模式（能量/角动量守恒）与近地 J2 长期摄动模式。
* **可见窗口**：给定地面站、最低仰角、时间区间；边界二分精化，
  收紧精度单调收敛、无跳变/丢失/乱序。
* **地影**：圆柱与锥形本影模型，输出进影/出影事件与累计光照时长。
* **批量并发**：可配置并发上限、队列深度、等待超时；结果确定性有序。
* **非物理输入显式失败**：e 越界、半长轴过低、迭代不收敛均返回可区分错误。
* **确定性与持久化**：相同输入+配置结果一致；场景落盘，重启后可重算。

## 文档

* [docs/conventions.md](docs/conventions.md) — 时间系统、坐标框架、单位
* [docs/perturbation.md](docs/perturbation.md) — 二体与 J2 摄动模型
* [docs/accuracy.md](docs/accuracy.md) — 精度契约、容差、TLE 一致性
* [docs/shadow.md](docs/shadow.md) — 圆柱/锥形地影模型与差异
* [docs/api.md](docs/api.md) — HTTP 接口与错误码

## 目录结构

```
app/
  timeutil.py    时间解析、JD、GMST
  elements.py    经典根数、开普勒方程、根数↔状态、物理校验
  tle.py         TLE 解析/校验/转换
  propagator.py  二体与 J2 解析传播器
  reference.py   独立 J2 RK4 参考传播器（测试基准）
  frames.py      J2000↔ECEF、WGS-84、ENU
  station.py     仰角与可见窗口（二分边界）
  shadow.py      太阳位置、圆柱/锥形地影、光照事件
  scheduler.py   有界并发队列
  storage.py     原子持久化
  service.py     业务编排
  main.py        FastAPI 路由
tests/           自动化测试
docs/            设计文档
```

## 已知限制

J2 解析模型为一阶平均根数模型（典型 LEO 72h 对独立数值参考 ≤25 km），
非 SGP4；忽略岁差章动极移、UT1−UTC、阻力/光压/三体摄动。详见
[docs/accuracy.md](docs/accuracy.md) 与 [docs/conventions.md](docs/conventions.md)。
