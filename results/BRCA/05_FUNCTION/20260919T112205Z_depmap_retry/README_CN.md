# DepMap访问重试

本轮问题：重试官方数据入口。

实际结果：本机DepMap网页、清单接口和Figshare历史接口均403；server165 requests也403。server165 curl得到HTTP200，但正文为DepMap Verification，含Cloudflare Turnstile人机验证组件，不是CSV数据。浏览器控制超时。

解释：确认至少该访问路径被人机验证阻挡，不能将HTTP200视为数据获取成功；其他403的精确原因未确定。

当前决定：ACCESS_BLOCKED，没有下载矩阵，没有计算依赖分数，不修改候选结论。

下一步：用户在普通浏览器完成网站正常验证后提供官方下载链接或可用文件，再在server165读取。

来源：https://depmap.org/portal/api/download/files
