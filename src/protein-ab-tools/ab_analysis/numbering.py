from typing import Literal, Optional
from anarci import anarci


def numbering(seq: str,
              name: Optional[str] = None,
              scheme: str = 'imgt',
              chain: Literal['H', 'L'] = 'H'):
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
