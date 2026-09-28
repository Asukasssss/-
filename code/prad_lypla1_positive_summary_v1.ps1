param([Parameter(Mandatory=$true)][string]$InputSummary,[Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$rows=Import-Csv -LiteralPath $InputSummary -Delimiter "`t"
$result=foreach($row in $rows){
 $n=[int]$row.n_cells
 $fraction=[double]$row.positive_fraction
 if($fraction -le 0 -or $fraction -gt 1){throw 'Invalid positive fraction'}
 $positive=$n*$fraction
 if([math]::Abs($positive-[math]::Round($positive)) -gt 1e-8){throw 'Positive count is not recoverable as integer'}
 [PSCustomObject]@{celltype=$row.celltype;n_all_cells=$n;n_positive_cells=[int][math]::Round($positive);positive_fraction=$fraction;positive_mean_log1p_CP10K=([double]$row.pooled_mean/$fraction);weighting='equal_cell';p_value='NA';q_value='NA';test_status='NOT_RUN'}
}
$result | Export-Csv -LiteralPath (Join-Path $OutputDir 'positive_cell_summary.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
$utf8=[Text.UTF8Encoding]::new($false)
$spec=@{version='prad_lypla1_positive_summary_v1';gene='LYPLA1';scope='PRAD cancer tissue, author normal versus malignant epithelium';formula='positive_mean_log1p_CP10K = pooled_mean_log1p_CP10K / positive_fraction';justification='Zero-count cells have log1p(CP10K)=0';weighting='equal cell, not patient equal';source='results/PRAD/06_EXTERNAL/20260928T031839Z_lypla1_umap_v2/epithelial_summary.tsv';statistics='NOT_RUN: no P/q, distribution test, or equivalence test';limitations='Conditional means only; not mean raw counts or mean unlogged CP10K; does not establish equal distributions';powershell=$PSVersionTable.PSVersion.ToString()}
[IO.File]::WriteAllText((Join-Path $OutputDir 'analysis_spec.json'),($spec|ConvertTo-Json -Depth 5),$utf8)
$validation=@{status='DONE';source_rows=$rows.Count;positive_counts_integral=$true;positive_fractions_valid=$true;individual_data_accessed=$false;malignant_minus_normal_mean=($result[1].positive_mean_log1p_CP10K-$result[0].positive_mean_log1p_CP10K)}
if($result[0].celltype -ne 'Normal epithelium' -or $result[1].celltype -ne 'Malignant epithelium'){throw 'Unexpected group ordering'}
[IO.File]::WriteAllText((Join-Path $OutputDir 'validation.json'),($validation|ConvertTo-Json),$utf8)
@([PSCustomObject]@{source=$spec.source;sha256=(Get-FileHash -LiteralPath $InputSummary -Algorithm SHA256).Hash.ToLower()},[PSCustomObject]@{source='code/prad_lypla1_positive_summary_v1.ps1';sha256=(Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLower()}) | Export-Csv -LiteralPath (Join-Path $OutputDir 'source_manifest.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
$result|Format-Table
