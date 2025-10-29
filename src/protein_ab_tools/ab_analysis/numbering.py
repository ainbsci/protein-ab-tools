from typing import Literal, Optional, List
from anarci import anarci


def run_numbering(
    seq: str,
    name: Optional[str] = None,
    scheme: str = 'imgt',
    chain: Literal['H', 'L'] = 'H',
    germline: bool = False,
    species: Optional[List[str]] = None,
):
    """
    Numbering an antibody sequence.
    """
    seq = seq.strip().replace('-', '')
    if name is None:
        name = f'{chain}-{scheme}'
    prep_seq = (name, seq)
    if chain == 'H':
        chain = ['H']
    elif chain == 'L':
        chain = ['K', 'L']
    result = anarci(
        [prep_seq],
        scheme=scheme,
        allow=chain,
        assign_germline=germline,
        allowed_species=species if species is not None else ['human', 'mouse'])
    if result[0][0] is None:
        raise ValueError(f"Invalid sequence: {seq}")
    return result


def extract_species(seq: str,
                    scheme: str = 'imgt',
                    chain: Literal['H', 'L'] = 'H'):
    result = run_numbering(seq, scheme=scheme, chain=chain)
    return result[0][0][1]['species']


def get_numbered_seq(seq: str,
                     scheme: str = 'imgt',
                     chain: Literal['H', 'L'] = 'H'):
    result = run_numbering(seq, scheme=scheme, chain=chain)
    return ''.join([aa for _, aa in result[0][0][0][0]])


def extract_regions(seq: str,
                    scheme: str = 'imgt',
                    chain: Literal['H', 'L'] = 'H'):
    """
    Extract regions from an antibody sequence.
    """
    result = run_numbering(seq, scheme=scheme, chain=chain)
    if scheme.lower() == 'imgt':
        breakpoint = {
            'fwr1': [1, 26],
            'cdr1': [27, 38],
            'fwr2': [39, 55],
            'cdr2': [56, 65],
            'fwr3': [66, 104],
            'cdr3': [105, 117],
            'fwr4': [118, 128]
        }
    elif scheme.lower() == 'aho':
        # https://plueckthun.bioc.uzh.ch/antibody/Numbering/NumFrame.html
        breakpoint = {
            'fwr1': [1, 26],
            'cdr1': [27, 40],
            'fwr2': [41, 57],
            'cdr2': [58, 68],
            'fwr3': [69, 106],
            'cdr3': [107, 138],
            'fwr4': [139, 149] if chain == 'H' else [139, 148]
        }
    elif scheme.lower() == 'chothia':
        breakpoint = {
            'fwr1': [1, 25] if chain == 'H' else [1, 23],
            'cdr1': [26, 32] if chain == 'H' else [24, 34],
            'fwr2': [33, 49] if chain == 'H' else [35, 49],
            'cdr2': [52, 56] if chain == 'H' else [50, 56],
            'fwr3': [63, 94] if chain == 'H' else [57, 88],
            'cdr3': [95, 102] if chain == 'H' else [89, 97],
            'fwr4': [103, 146] if chain == 'H' else [98, 146]
        }
    elif scheme.lower() == 'kabat':
        breakpoint = {
            'fwr1': [1, 30] if chain == 'H' else [1, 23],
            'cdr1': [31, 35] if chain == 'H' else [24, 34],
            'fwr2': [36, 49] if chain == 'H' else [35, 49],
            'cdr2': [50, 65] if chain == 'H' else [50, 56],
            'fwr3': [66, 94] if chain == 'H' else [57, 88],
            'cdr3': [95, 102] if chain == 'H' else [89, 97],
            'fwr4': [103, 146] if chain == 'H' else [98, 146]
        }
    else:
        raise ValueError('Invalid numbering scheme')
    chain = 'vh' if chain == 'H' else 'vl'
    numbered_seq = result[0][0][0][0]
    regions = {f'{chain}_{k}': '' for k in breakpoint.keys()}
    for (number, _), aa in numbered_seq:
        region = next(
            (k for k, v in breakpoint.items() if v[0] <= number <= v[1]), None)
        if region is not None:
            regions[f'{chain}_{region}'] += aa
    return regions


if __name__ == '__main__':

    right_h = {
        'vh_fwr1': 'QVQLVESGG-GVVQPGRSLRLDCKAS',
        'vh_cdr1': 'GITF----SNSG',
        'vh_fwr2': 'MHWVRQAPGKGLEWVAV',
        'vh_cdr2': 'IWYD--GSKR',
        'vh_fwr3': 'YYADSVK-GRFTISR-NSKNTLFLQMNSLRAEDTAVYYC',
        'vh_cdr3': 'ATN-------DDY',
        'vh_fwr4': 'WGQGTLVTTVSS'
    }
    right_l = {
        'vl_fwr1': 'EIVLTQSPATLSLSPGERATLSCRAS',
        'vl_cdr1': 'QSV------SGY',
        'vl_fwr2': 'LAWYQQKPGQAPRLLIY',
        'vl_cdr2': 'DA-------S',
        'vl_fwr3': 'NRATGIP-ARFSGSG--SGTDFTLTISSLEPEDFAVYYC',
        'vl_cdr3': 'QQSSN----WPRT',
        'vl_fwr4': 'FGQGTKVEIK'
    }
    expected_h = extract_regions(
        'QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTTVSS',
        chain='H')
    assert expected_h == right_h
    expected_l = extract_regions(
        'EIVLTQSPATLSLSPGERATLSCRASQSVSGYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQSSNWPRTFGQGTKVEIK',
        chain='L')
    assert expected_l == right_l
