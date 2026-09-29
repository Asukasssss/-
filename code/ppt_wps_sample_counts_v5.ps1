$ErrorActionPreference='Stop'
$repo=Split-Path $PSScriptRoot -Parent
$out=Join-Path $repo 'results/BRCA/07_INTEGRATION/20260929T_wps_counts_v5'
New-Item -ItemType Directory -Force -Path $out | Out-Null
$view='D:/CodexData/visualizations/2026/09/29/01a0ec94-4db2-7ec1-8e34-4f52c52a8ce1/CAMP_case_art_v4'
$target=Join-Path $view 'CAMP_完整组会汇报_生化插图精修_25页.pptx'
$app=New-Object -ComObject KWPP.Application
$deck=$null
for($i=1;$i -le $app.Presentations.Count;$i++) { $p=$app.Presentations.Item($i);if($p.FullName -eq $target.Replace('/','\')){$deck=$p;break} }
if($null -eq $deck){throw 'Expected open WPS deck not found; no documents changed'}
$slide=$deck.Slides.Item(4)
$updates=@(
 @{id=12;text="108 份标本`r其中配对：45对 / 90份";top=267;size=11},
 @{id=17;text="66 份标本`r其中配对：33对 / 66份";top=267;size=11},
 @{id=22;text="39 份标本`r其中配对：11对 / 22份";top=267;size=11},
 @{id=27;text="137 份标本`r其中配对：43对 / 86份";top=400.2;size=10.5},
 @{id=32;text="队列3：114份，配对17位`r队列4：71份，配对12位";top=400.2;size=10.5},
 @{id=37;text="80 份：74 肿瘤 + 6 参照`r0 对 · 跨来源非配对";top=400.2;size=10.5}
)
foreach($u in $updates){
 $shape=$slide.Shapes.Item($u.id);$shape.Top=[single]$u.top;$shape.Height=[single]34;$shape.Width=[single]148
 $shape.TextFrame.TextRange.Text=$u.text
 $shape.TextFrame.TextRange.Font.Name='Microsoft YaHei';$shape.TextFrame.TextRange.Font.NameFarEast='微软雅黑';$shape.TextFrame.TextRange.Font.Size=[single]$u.size
 $shape.TextFrame.TextRange.ParagraphFormat.SpaceBefore=0;$shape.TextFrame.TextRange.ParagraphFormat.SpaceAfter=2
}
$slide.Shapes.Item(6).TextFrame.TextRange.Text='当前分析范围：标本份数与同患者肿瘤—参照配对数分别列示'
$foot=$slide.Shapes.Item(39)
$foot.Top=501;$foot.Height=36;$foot.Width=810
$foot.TextFrame.TextRange.Text="配对按患者计；ccRCC为多区域标本，17 / 12对分别对应17 / 12位患者，两队列不合并。`rBRCA含1份组织标签冲突；PDAC 39份中11对可配对。各具体分析按实际可用数计。"
$foot.TextFrame.TextRange.Font.Name='Microsoft YaHei';$foot.TextFrame.TextRange.Font.NameFarEast='微软雅黑';$foot.TextFrame.TextRange.Font.Size=9
$foot.TextFrame.TextRange.ParagraphFormat.SpaceBefore=0;$foot.TextFrame.TextRange.ParagraphFormat.SpaceAfter=1
$deck.Save()
$slide.Export((Join-Path $view 'rendered/幻灯片4.PNG'),'PNG',1600,900)
$deck.SaveAs((Join-Path $out 'CAMP_样本数补充_25页预览.pdf'),32)
# SaveAs PDF can change the active filename on some WPS versions; retain the editable PPT.
$deck.SaveAs($target,24)
Copy-Item -LiteralPath $target -Destination (Join-Path $out 'CAMP_样本数补充_25页.pptx') -Force
Copy-Item -LiteralPath (Join-Path $view 'rendered/幻灯片4.PNG') -Destination (Join-Path $out '04_队列样本与配对.png') -Force
Copy-Item -LiteralPath (Join-Path $out 'CAMP_样本数补充_25页预览.pdf') -Destination (Join-Path $view 'CAMP_生化插图精修_预览.pdf') -Force
$app.ActiveWindow.View.GotoSlide(4)
Write-Output 'Updated slide 4 directly in WPS; saved PPT and refreshed preview.'
