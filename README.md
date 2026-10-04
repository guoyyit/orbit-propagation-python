# Orbit Propagation Service (Python)

纯本地的卫星轨道预报与地面站可见性计算服务。无任何外部服务依赖，输入确定
则结果可逐位复现。

## 功能

- TLE 与经典轨道六根数两种输入，字段越界、校验和、格式错误返回可区分错误
  码与出错位置。
- 任意时刻位置/速度外推；时间支持 ISO 8601 与 Unix 秒。近地轨道计入 J2
  摄动（定步长 RK4）；无摄动模式能量与角动量守恒。
- 给定地面站、最低仰角与时间区间的可见窗口序列；边界精度可配置，收紧后
  单调收敛、不跳变/丢失/乱序。
- 圆柱与锥形两种地影模型，输出进/出影事件与累计光照时长。
- 批量并发查询，受并发上限与队列深度、等待超时约束；结果保序且确定。
- 非物理输入（偏心率越界、半长轴过低/穿地、迭代不收敛）显式失败。
- 场景以原子 JSON 持久化，重启后可重算且结果一致。

## 快速开始

```bash
uv sync --group dev
uv run uvicorn main:app --reload
uv run pytest
```

```bash
curl -X PUT localhost:8000/orbits/iss -H 'content-type: application/json' -d '{
  "format":"coe","angles":"degrees",
  "semi_major_axis":6794900,"eccentricity":0.0015,
  "inclination":51.64,"raan":120.0,
  "argument_of_periapsis":60.0,"mean_anomaly":0.0,
  "epoch":"2023-11-14T22:13:20Z"}'
curl -X POST localhost:8000/propagate -H 'content-type: application/json' -d '{
  "orbit":"iss","time":"2023-11-15T00:00:00Z"}'
```

## 文档

- [时间与坐标约定](docs/conventions.md)
- [摄动模型与数值方法](docs/perturbation.md)
- [精度约定与容差](docs/accuracy.md)
- [地影模型](docs/eclipse.md)
- [HTTP API 与错误码](docs/api.md)
- [已知限制](docs/limitations.md)
