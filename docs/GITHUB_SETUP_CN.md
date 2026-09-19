# GitHub 发布与第二台电脑

用户已指定公共仓库 https://github.com/Asukasssss/- 。本机 Git 使用 Asukasssss 身份认证；这不代表 Codex GitHub 插件已切换账号。

当前电脑在仓库目录执行（身份认证由用户正常登录GitHub完成，不能把token写在URL/文件）：

```powershell
git remote add origin https://github.com/Asukasssss/-.git
git push -u origin main
```

若origin已存在，先git remote -v核对，不能盲目替换。远程若已有提交，先读取并整合，不能force push。

另一台电脑：

```powershell
git clone https://github.com/Asukasssss/-.git camp-cancer-collaboration
cd camp-cancer-collaboration
```

在Codex里打开此目录，把START_ACCOUNT_B_CN.md交给它。使用有仓库权限的GitHub身份；Codex账号与GitHub仓库权限分别配置，不需要互相分享密码。未获权限时由用户完成授权，不在交接包存凭据。

若GitHub还没就绪，可以直接传递已生成的本地zip或git bundle。只有DEPLOYMENT_STATUS_CN.md确认服务器副本就绪后，才从server165下载。bundle可用 `git clone <本地bundle路径> camp-cancer-collaboration`；之后添加真实origin。此临时方案不代表GitHub已上线。
