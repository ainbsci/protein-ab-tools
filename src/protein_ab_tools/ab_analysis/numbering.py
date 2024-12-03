from typing import Literal, Optional
from anarci import anarci


def run_numbering(
    seq: str,
    name: Optional[str] = None,
    scheme: str = 'imgt',
    chain: Literal['H', 'L'] = 'H',
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
    result = anarci([prep_seq], scheme=scheme, allow=chain)
    if result[0][0] is None:
        raise ValueError(f"Invalid sequence: {seq}")
    return result


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
            'fwr1': [1, 27],
            'cdr1': [27, 39],
            'fwr2': [39, 56],
            'cdr2': [56, 66],
            'fwr3': [66, 105],
            'cdr3': [105, 118],
            'fwr4': [118, 129]
        }
    elif scheme.lower() == 'aho':
        # https://plueckthun.bioc.uzh.ch/antibody/Numbering/NumFrame.html
        breakpoint = {
            'fwr1': [1, 27],
            'cdr1': [27, 41],
            'fwr2': [41, 58],
            'cdr2': [58, 69],
            'fwr3': [69, 107],
            'cdr3': [107, 139],
            'fwr4': [139, 150] if chain == 'H' else [139, 149]
        }
    else:
        raise ValueError('Invalid numbering scheme')
    chain = 'vh' if chain == 'H' else 'vl'
    numbered_seq = result[0][0][0][0]
    regions = {f'{chain}_{k}': '' for k in breakpoint.keys()}
    for (number, _), aa in numbered_seq:
        region = next(
            (k for k, v in breakpoint.items() if v[0] <= number < v[1]), None)
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
