# 病理图取数补查：16个切片详情及关联仓库

日期：2026-09-28。接续前批取数结果，不新增表达分析或叠加。

## 实际补查

|入口|读取结果|结论边界|
|---|---|---|
|STTS0000650–STTS0000665，共16个切片详情页|全部成功读取服务器返回的页面状态；均记H&E Staining，未发现除站点标识外的img标签，详情数据未给出图像地址|页面无图片链接不证明原图不存在；不是只查第一张|
|StereoMM仓库分支与release|只有StereoMM分支，未列出release；此前递归文件树没有配套CRC扫描图|示例image.tif路径不是可下载实物|
|StereoMM issue #1及作者回复|提问者请求配套数据；作者指向STOmics官网小鼠肾脏H&E示例|没有提供当前CRC切片；不拿小鼠数据替代|
|原文Peer Review File（MOESM2）全文文本链接及相关段落|数据说明仍指向GSA与STT0000036，分析绘图代码可向作者索取；未见额外病理原图地址|未对审稿文件全部插图作病理复核|
|Additional Supplementary Information说明（MOESM3）|四份数据表分别为临床信息、单细胞差异基因、注释基因、配体受体结果|不是一个遗漏的病理图像附件目录|
|CNP0002432目录|Single_Cell、CRC_SingleCell.rds、CRC_Stereo.rds；Stereo RDS约1.04GB|属于另一研究的CRC数据；未下载或检查RDS内部，不能声称已排除其图像内容，也不能当成当前16切片|

## 直接来源

- [16切片列表](https://db.cngb.org/stomics/project/STT0000036/spatialGeneExpression)
- [第一张切片详情](https://db.cngb.org/stomics/tissue_section/STTS0000650)，其余链接保存在tissue_page_followup.json。
- [StereoMM作者答复](https://github.com/STOmics/StereoMMv1/issues/1#issuecomment-2510409465)
- [StereoMM仓库](https://github.com/STOmics/StereoMMv1/tree/StereoMM)
- [原论文补充材料](https://www.nature.com/articles/s41467-024-54710-3#Sec25)
- [另一CRC研究CNP0002432文件目录](https://ftp.cngb.org/pub/CNSA/data4/CNP0002432/)

## 本轮结论

补充图2a的H&E图已经取得，但此次扩大公开入口核查仍未取得当前STT0000036的原始高分辨率切片与配准文件。保持PARTIAL，不把“未找到”升级成“作者未公开”或“不存在”。若要继续这批的精确叠加，仍需明确图像下载位置、图像—表达标本映射及变换参数。未发送邮件或GitHub留言。

## 复现与校验

code/coad/search_stereoseq_histology_followup.py从仓库根目录运行，读取公开切片详情、关联FTP目录、StereoMM分支/release及论文XML；结果暂存runtime，不下载矩阵。tissue_page_followup.json保留本轮16页公开处理信息与读取状态，不含患者表达值。额外审稿文件与DOCX在server165解包并只读检查，源文件不进仓库。所有16页成功返回H&E元数据；网页核查不等于全文件存储审计。
