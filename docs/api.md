# HTTP API

本地服务，默认 `uvicorn app.main:app`，无任何外部网络依赖。

## 端点

### `GET /health`
健康检查，返回 `{"status":"ok"}`。

### `POST /api/v1/propagate`
对一组时刻外推位置/速度（J2000）。

```json
{
  "orbit": {
    "orbit_id": "o1",
    "perturbation": "j2",
    "elements": {
      "epoch": "2024-03-20T00:00:00Z",
      "semi_major_axis_m": 7000000.0,
      "eccentricity": 0.001,
      "inclination_rad": 0.9013,
      "raan_rad": 4.318,
      "arg_perigee_rad": 0.698,
      "mean_anomaly_rad": 0.1745
    }
  },
  "times": ["2024-03-20T00:30:00Z", 1710894600]
}
```

也可用 `tle` 替代 `elements`：

```json
"tle": {"line0": "ISS", "line1": "1 ...", "line2": "2 ..."}
```

`times` 元素既可是 ISO 8601 字符串，也可是 Unix 秒。`perturbation` 取
`none` 或 `j2`（默认 `j2`）。

### `POST /api/v1/visibility`
地面站可见窗口。

```json
{
  "orbit": { "...": "同 propagate" },
  "station": {"name": "gs", "latitude_deg": 0.0, "longitude_deg": 0.0, "altitude_m": 0.0},
  "start": "2024-03-20T00:00:00Z",
  "end":   "2024-03-20T12:00:00Z",
  "min_elevation_deg": 5.0,
  "boundary_tol_s": 1.0,
  "sample_step_s": 60.0
}
```

返回按时间排序的窗口数组（`start`/`end`/`max_elevation_deg`）。

### `POST /api/v1/illumination`
地影事件与光照时长。`model` 取 `cylindrical` 或 `conical`（默认）。
返回 `events`（`entering_shadow`/`leaving_shadow`）、
`illuminated_seconds`、`shadowed_seconds`。

### `POST /api/v1/batch`
一次提交多个轨道（与多个地面站），在配置的并发上限内执行、超出排队。
结果按输入的确定性顺序返回。并发参数通过服务配置
（`SchedulerConfig`: `max_concurrency` / `max_queue_depth` / `wait_timeout_s`）。

### `POST /api/v1/scenarios` / `GET /api/v1/scenarios/{id}`
注册/读取场景。保存为本地 JSON（原子写+fsync），服务重启后仍可读取，
重算结果与重启前一致（纯函数计算，无随机/时钟依赖）。

## 错误码

所有错误返回 `{code, message, location}`，HTTP 状态：

| code | HTTP | 含义 |
|------|------|------|
| `VALIDATION_ERROR` | 400 | 通用请求校验错误 |
| `TLE_FORMAT_ERROR` | 400 | TLE 行结构/长度/非法字符 |
| `TLE_CHECKSUM_ERROR` | 400 | NORAD 模 10 校验和不符 |
| `TLE_FIELD_ERROR` | 400 | TLE 字段越界 |
| `FIELD_OUT_OF_RANGE` | 400 | 数值字段越界（如 e≥1） |
| `NON_PHYSICAL_ORBIT` | 400 | 近地点低于地表等非物理轨道 |
| `TIME_FORMAT_ERROR` | 400 | 时间字符串/时间戳非法 |
| `ITERATION_NOT_CONVERGED` | 400 | 开普勒方程等迭代不收敛 |
| `QUEUE_FULL` | 503 | 有界队列已满 |
| `QUEUE_TIMEOUT` | 503 | 排队等待超时 |
| `NOT_FOUND` | 404 | 场景不存在 |
| `STORAGE_ERROR` | 500 | 持久化读写失败 |

## 精度与边界参数

* `boundary_tol_s`：窗口/地影边界二分精度（秒），收紧单调收敛。
* `sample_step_s`：扫描步长（秒），决定最短可检测事件，默认 60。

模型精度、时间/坐标约定、地影差异见同目录其它文档。
