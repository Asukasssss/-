$ErrorActionPreference = 'Stop'
$root = Join-Path $PSScriptRoot '../results/BRCA/07_INTEGRATION/20260929T_illustrated_main_v2'
$root = (Resolve-Path -LiteralPath $root).Path
$render = Join-Path $root 'rendered'
New-Item -ItemType Directory -Force -Path $render | Out-Null
$app = New-Object -ComObject KWPP.Application
$deck = $app.Presentations.Open((Join-Path $root 'CAMP_泛癌与LYPLA1_插图精修可编辑版.pptx'),0,0,0)
try {
 if ($deck.Slides.Count -ne 18) { throw 'Unexpected slide count' }
 $deck.Export($render,'PNG',1600,900)
 $deck.SaveAs((Join-Path $root 'CAMP_插图精修版_预览.pdf'),32)
} finally { $deck.Close() }
Write-Output 'WPS_EXPORT_DONE'
