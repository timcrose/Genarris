from evaluation.fhi_aims_utils import get_fhi_aims_extractor

def fhi_aims_extract(inst):
    sname = "fhi_aims_extract"
    calculation_dir = inst.get(sname, "calculation_dir")
    fhi_aims_extractor = get_fhi_aims_extractor(inst, sname)
    return fhi_aims_extractor.extract_batch(calculation_dir)
