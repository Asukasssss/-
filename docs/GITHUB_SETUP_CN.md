# GitHub 发布与第二台电脑

本地仓库已准备；若远程尚未存在，先在DengLeiO名下创建私有仓库camp-cancer-collaboration，不初始化README。不要使用已有无关的公共仓库。

当前电脑在仓库目录执行（身份认证由用户正常登录GitHub完成，不能把token写在URL/文件）：

```powershell
git remote add origin https://github.com/DengLeiO/camp-cancer-collaboration.git
git push -u origin main
```

若origin已存在，先git remote -v核对，不能盲目替换。远程若已有提交，先读取并整合，不能force push。

另一台电脑：

```powershell
git clone https://github.com/DengLeiO/camp-cancer-collaboration.git
cd camp-cancer-collaboration
```

在Codex里打开此目录，把START_ACCOUNT_B_CN.md交给它。使用有仓库权限的GitHub身份；Codex账号与GitHub仓库权限分别配置，不需要互相分享密码。未获权限时由用户完成授权，不在交接包存凭据。

若GitHub还没就绪，可从server165交接目录下载zip或git bundle，先读取上下文。bundle可用 `git clone <本地bundle路径> camp-cancer-collaboration`；之后添加真实origin。此临时方案不代表GitHub已上线。
