$ErrorActionPreference='Stop'
$repo=Split-Path $PSScriptRoot -Parent
$o=Join-Path $repo 'results/BRCA/07_INTEGRATION/20260929T_reaction_inplace_v8';New-Item -ItemType Directory -Force $o | Out-Null
$v='D:/CodexData/visualizations/2026/09/29/01a0ec94-4db2-7ec1-8e34-4f52c52a8ce1/CAMP_reaction_inplace_v8';New-Item -ItemType Directory -Force $v | Out-Null
$a=New-Object -ComObject KWPP.Application;$d=$a.ActivePresentation
if($d.Slides.Count -ne 26){throw 'Expected 26-slide deck'}
$s=$d.Slides.Item(9);$src=$d.Slides.Item(10)
foreach($id in @(51,50,49,48,47,46,45)){$s.Shapes.Item($id).Delete()}
$s.Shapes.Item(6).Height=[single]102
function T($x,$y,$w,$h,$text,$size,$bold){
 $q=$s.Shapes.AddTextbox(1,[single]$x,[single]$y,[single]$w,[single]$h)
 $q.TextFrame.MarginLeft=0;$q.TextFrame.MarginRight=0;$q.TextFrame.MarginTop=0;$q.TextFrame.MarginBottom=0
 $q.TextFrame.WordWrap=0;$r=$q.TextFrame.TextRange;$r.Text=$text;$r.Font.Name='Microsoft YaHei';$r.Font.NameFarEast='微软雅黑';$r.Font.Size=[single]$size;$r.Font.Color.RGB=3888384
 $r.ParagraphFormat.SpaceAfter=0
}
function Pic($id,$x,$y,$height){
 $src.Shapes.Item($id).Copy();$p=$s.Shapes.Paste().Item(1);$p.LockAspectRatio=-1;$p.Height=[single]$height;$p.Left=[single]$x;$p.Top=[single]$y
}
foreach($row in @(@{y=193;ids=@(9,10,11);left='LPC(16:0)';right='GPC + 棕榈酸';role='水解 LPC，释放棕榈酸'},@{y=241;ids=@(20,21,22);left='含油酰基 PC';right='2-酰基 LPC + 油酸';role='切下 sn-1 油酰基，释放油酸'})){
 $y=$row.y;Pic $row.ids[0] 516 ($y+2) 29;Pic $row.ids[1] 651 ($y+1) 20;Pic $row.ids[2] 757 ($y+1) 30
 T 544 ($y+9) 100 16 $row.left 11 -1
 T 677 ($y+1) 74 14 'LYPLA1' 10 -1
 T 790 ($y+10) 120 16 $row.right 10.3 -1
 $l=$s.Shapes.AddLine(643,[single]($y+25),746,[single]($y+25));$l.Line.ForeColor.RGB=3888384;$l.Line.Weight=1.2;$l.Line.EndArrowheadStyle=3
 T 545 ($y+33) 355 12 $row.role 8.5 0
}
$d.Slides.Item(10).Delete()
for($n=10;$n -le $d.Slides.Count;$n++){$sl=$d.Slides.Item($n);for($j=1;$j -le $sl.Shapes.Count;$j++){$q=$sl.Shapes.Item($j);if($q.Left -gt 850 -and $q.Top -gt 490 -and $q.HasTextFrame){if($q.TextFrame.TextRange.Text -match '^\d+$'){$q.TextFrame.TextRange.Text=[string]$n}}}}
$ppt=Join-Path $v 'CAMP_反应示意原位替换_25页.pptx';$d.SaveAs($ppt,24)
$s.Export((Join-Path $v '09_原位反应示意.png'),'PNG',1600,900)
$d.SaveAs((Join-Path $v 'CAMP_25页预览.pdf'),32);$d.SaveAs($ppt,24)
Get-ChildItem -LiteralPath $v -File | Where-Object {$_.Name -notlike '~*'} | Copy-Item -Destination $o -Force
$a.ActiveWindow.View.GotoSlide(9)
Write-Output 'Replaced slide 9 reaction panel, removed standalone page, saved 25 slides.'
