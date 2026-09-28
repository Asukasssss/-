"""Conservative author-label identity resolution; never infer malignancy from tissue."""
MALIGNANT='Malignant (author)'
NORMAL='Normal epithelial (author)'
UNRESOLVED='Ductal (unresolved)'

def resolve_identity(cohort,broad,fine=None):
    fine='' if fine is None else str(fine).strip()
    if cohort=='GSE242230':
        if fine=='Normal Epithelial':return NORMAL,'AUTHOR_NORMAL','Specific author label overrides broad Malignant label'
        if fine in ['Malignant - Classical','Malignant - Basal']:return MALIGNANT,'AUTHOR_MALIGNANT','Explicit author malignant subtype;not independent CNV validation'
        if broad=='Malignant':return UNRESOLVED,'UNRESOLVED','Broad Malignant without a recognized specific identity'
    if broad in ['Ductal','Ductal cell']:return UNRESOLVED,'UNRESOLVED','GEO supplied broad ductal label only;no per-cell malignant call'
    return broad,'OTHER_LINEAGE','Original nonductal lineage retained'
