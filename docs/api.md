# HTTP API

所有接口接收/返回 JSON。错误统一形如：

```json
{"error": {"code": "non_physical_orbit", "message": "...", "location": "eccentricity"}}
```

## 错误码

| code | HTTP | 含义 |
| --- | --- | --- |
| `tle_format_error` | 422 | TLE 结构/列宽/数字格式错误，`location` 指向行列 |
| `tle_checksum_error` | 422 | TLE 模 10 校验和错误，`location=lineN:69` |
| `field_range_error` | 422 | 字段越界（角度、经纬度、容差等） |
| `time_format_error` | 422 | 时间既非 ISO 8601 也非合法 Unix 秒 |
| `non_physical_orbit` | 422 | 偏心率越界、近地点穿地、非束缚轨道等 |
| `convergence_error` | 422 | 开普勒/求根/积分迭代不收敛或出现非有限值 |
| `queue_overflow` | 503 | 超出并发+队列上限 |
| `wait_timeout` | 504 | 排队等待超过 `wait_timeout_seconds` |
| `scene_not_found` | 404 | 引用了未注册的轨道/地面站 |
| `scene_conflict` | 409 | 注册重名场景且未用 PUT 覆盖 |

## 端点

- `GET /health`：健康检查与已注册场景列表。
- `POST /orbits/{name}` / `PUT /orbits/{name}`：注册（POST 禁重名，PUT 覆盖）。
- `POST /stations/{name}` / `PUT /stations/{name}`：注册地面站。
- `GET /scenes` / `DELETE /scenes/{name}`：列举/删除场景。
- `POST /propagate`：体 `{"orbit": <轨道或引用>, "time": ..., "propagator": {...}}`。
- `POST /visibility`：体含 `orbit`、`station`、`start`、`end`、
  `minimum_elevation_deg`、`boundary_tolerance_seconds`。
- `POST /illumination`：体含 `orbit`、`start`、`end`、
  `shadow_model`（`cylindrical`/`conical`）、`boundary_tolerance_seconds`。
- `POST /batch`：体 `{"queries": [...]}`，每项含 `operation`
  （`propagate`/`visibility`/`illumination`）及对应参数。

轨道对象可以是内联定义，也可以是已注册名字的字符串，或 `{"ref": "name"}`。

## 并发与队列

`SchedulerConfig`：`max_concurrency`（默认 4）、`max_queue_depth`（默认 32）、
`wait_timeout_seconds`（默认 30）。批量结果严格按提交顺序返回，与实际完成
先后无关；同一输入在任意并发下逐位一致。

## 持久化

注册场景以原子写入的 JSON 保存在 `data/scenes.json`（可配置）。重启后场景
仍可用，且因存储的只是输入、重算路径确定，重启前后结果一致。

## 本地运行

```bash
uv sync --group dev
uv run uvicorn main:app --host 127.0.0.1 --port 8000
uv run pytest
```
