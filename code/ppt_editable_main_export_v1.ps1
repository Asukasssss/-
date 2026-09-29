$ErrorActionPreference = 'Stop'
$root = Join-Path $PSScriptRoot '../results/BRCA/07_INTEGRATION/20260929T_editable_main_v1'
$root = (Resolve-Path -LiteralPath $root).Path
$render = Join-Path $root 'rendered'
New-Item -ItemType Directory -Force -Path $render | Out-Null
$app = New-Object -ComObject KWPP.Application
$deck = $app.Presentations.Open((Join-Path $root 'CAMP_泛癌发现与LYPLA1_可编辑汇报.pptx'),0,0,0)
try {
 if ($deck.Slides.Count -ne 18) { throw 'Unexpected slide count' }
 $deck.Export($render,'PNG',1600,900)
 $deck.SaveAs((Join-Path $root 'CAMP_泛癌发现与LYPLA1_预览.pdf'),32)
} finally { $deck.Close() }
Write-Output 'WPS_EXPORT_DONE'
