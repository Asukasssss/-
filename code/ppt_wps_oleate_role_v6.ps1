$ErrorActionPreference='Stop'
$repo=Split-Path $PSScriptRoot -Parent
$out=Join-Path $repo 'results/BRCA/07_INTEGRATION/20260929T_wps_oleate_v6'
New-Item -ItemType Directory -Force -Path $out | Out-Null
$view='D:/CodexData/visualizations/2026/09/29/01a0ec94-4db2-7ec1-8e34-4f52c52a8ce1/CAMP_case_art_v4'
$target=Join-Path $view 'CAMP_完整组会汇报_生化插图精修_25页.pptx'
$app=New-Object -ComObject KWPP.Application
$d=$app.ActivePresentation
if($d.FullName -ne $target.Replace('/','\')){throw 'Current WPS document differs; no change made'}
$s=$d.Slides.Item(9)
$q=$s.Shapes.Item(50)
$q.Left=[single]524.16;$q.Top=[single]245;$q.Width=[single]370;$q.Height=[single]17
$q.TextFrame.TextRange.Text='含油酰基 PC → 2-酰基 LPC + 油酸'
$q.TextFrame.TextRange.Font.Size=[single]11.2
$q.TextFrame.TextRange.Font.Name='Microsoft YaHei';$q.TextFrame.TextRange.Font.NameFarEast='微软雅黑'
$q=$s.Shapes.Item(51)
$q.Left=[single]524.16;$q.Top=[single]266;$q.Width=[single]370;$q.Height=[single]15
$q.TextFrame.TextRange.Text='油酸：sn-1 位水解释放的脂肪酸产物'
$q.TextFrame.TextRange.Font.Size=[single]10
$q.TextFrame.TextRange.Font.Name='Microsoft YaHei';$q.TextFrame.TextRange.Font.NameFarEast='微软雅黑'
$d.Save()
$s.Export((Join-Path $view 'rendered/幻灯片9.PNG'),'PNG',1600,900)
$d.SaveAs((Join-Path $out 'CAMP_油酸作用补充_预览.pdf'),32)
$d.SaveAs($target,24)
Copy-Item -LiteralPath $target -Destination (Join-Path $out 'CAMP_油酸作用补充_25页.pptx') -Force
Copy-Item -LiteralPath (Join-Path $view 'rendered/幻灯片9.PNG') -Destination (Join-Path $out '09_油酸反应说明.png') -Force
Copy-Item -LiteralPath (Join-Path $out 'CAMP_油酸作用补充_预览.pdf') -Destination (Join-Path $view 'CAMP_生化插图精修_预览.pdf') -Force
$app.ActiveWindow.View.GotoSlide(9)
Write-Output 'Saved WPS slide 9 with oleate reaction and product role.'
