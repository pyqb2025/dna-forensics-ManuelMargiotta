"""A Profiler provides methods to use a DNA string for profiling suspects."""


class Profiler:
    """Encapsulate a DNA sequence.

    >>> p = Profiler('CTAGATAGATAGATAGATGACTA')
    >>> p.longest_run('AGAT')
    4
    >>> p.match_suspect('Ada', {'AGAT': 4})
    True
    >>> p.match_suspect('Bia', {'AGAT': 5})
    False

    Using a DNA database (in CSV format):

    >>> with open('sequence.txt') as seq_file:
    ...     s = seq_file.read()
    ...     p = Profiler(s)
    ...     with open('data.csv') as dna_db:
    ...             keys = dna_db.readline().strip().split(',')[1:]
    ...             fpr = {k: 0 for k in keys}
    ...             results = []
    ...             for row in dna_db:
    ...                     cols = row.strip().split(',')
    ...                     for i, k in enumerate(keys):
    ...                             fpr[k] = int(cols[i+1])
    ...                     if p.match_suspect(cols[0], fpr):
    ...                             results.append(cols[0] + " guilty")
    ...                     else:
    ...                             results.append(cols[0] + " innocent")
    ...             results
    ['Andrew innocent', 'Athena innocent', 'Brian innocent', 'Chad innocent', 'David innocent', 'Doug innocent', 'Erin guilty', 'Ian innocent', 'Jelle innocent', 'Kareem innocent', 'Meredith innocent', 'Rodrigo innocent', 'Tara innocent', 'Teagan innocent', 'Valerie innocent']
    """

    def __init__(self, sequence: str):
        """Create a Profiler with a sequence.

        >>> p = Profiler('AAGCT')
        >>> p.seq
        'AAGCT'
        """
        self.seq = sequence

    def longest_run(self, subseq: str) -> int:
        """Return the longest number of repetitions of subseq in the encapsulated DNA sequence.

        >>> p = Profiler('AACCCTGCGCGCGCGCGATCTATCTATCTATCTATCCAGCATTAGCTAGCATCAAGATAGATAGATGAATTTCGAAATGAATGAATGAATGAATGAATGAATG')
        >>> p.longest_run('AGAT')
        3
        >>> p.longest_run('AATG')
        7
        >>> p.longest_run('TATC')
        4
        >>> p = Profiler('CCAGATAGATAGATAGATAGATAGATGTCACAGGGATGCTGAGGGCTGCTTCGTACGTACTCCTGATTTCGGGGATCGCTGACACTAATGCGTGCGAGCGGATCGATCTCTATCTATCTATCTATCTATCCTATAGCATAGACATCCAGATAGATAGATC')
        >>> p.longest_run('AGAT')
        6
        >>> p.longest_run('AATG')
        1
        >>> p.longest_run('TATC')
        5
        """
        c=0
        longest=0
        i=0
        while (i<len(self.seq)):
            if self.seq[i:i+len(subseq)]==subseq:
                c+=1
                i+=len(subseq)
                if c>longest:
                    longest=c
            else:    
                i+=1
                c=0
        return longest


    def match_suspect(self, suspect_name: str, dna_fpr: dict[str, int]) -> bool:
        """True if the dna_fpr associated to suspect_name can be found exactly in the DNA sequence. 

        >>> p = Profiler('AGACGGGTTACCATGACTATCTATCTATCTATCTATCTATCTATCTATCACGTACGTACGTATCGAGATAGATAGATAGATAGATCCTCGACTTCGATCGCAATGAATGCCAATAGACAAAA')
        >>> p.match_suspect('Cain', {'AGAT':5, 'AATG':2, 'TATC':8})
        True
        >>> p.match_suspect('Abel', {'AGAT':3, 'AATG':7, 'TATC':4})
        False
        """
        bool_list=[]
        suspect_dict:dict[str,int]={}
        for key in dna_fpr.keys():
            suspect_dict.update({key:self.longest_run(key)}) 
        for key in dna_fpr.keys():
            bool_list.append(suspect_dict[key]==dna_fpr[key])
        return all(bool_list)    

