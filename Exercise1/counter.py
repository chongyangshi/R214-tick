from collections import Counter
class TermCounts():
    '''
        TermCounts counts the term counts for each
    '''
    def __init__(self, go, annots):
        '''
            Initialise the counts and
        '''
        # Backup
        self._go = go

        # Initialise the counters
        self._counts = Counter()
        self._aspect_counts = Counter()

        # Fill the counters...
        self._count_terms(go, annots)

    def _count_terms(self, go, annots):
        '''
            Fills in the counts and overall aspect counts.
        '''
        for x in annots:
            # Extract term information
            go_id = annots[x]['GO_ID']
            namespace = go[go_id].namespace

            self._counts[go_id] += 1
            rec = go[go_id]
            parents = rec.get_all_parents()
            for p in parents:
                self._counts[p] += 1

            self._aspect_counts[namespace] += 1

    def get_count(self, go_id):
        '''
            Returns the count of that GO term observed in the annotations.
        '''
        return self._counts[go_id]

    def get_total_count(self, aspect):
        '''
            Gets the total count that's been precomputed.
        '''
        return self._aspect_counts[aspect]

    def get_term_freq(self, go_id):
        '''
            Returns the frequency at which a particular GO term has
            been observed in the annotations.
        '''
        try:
            namespace = self._go[go_id].namespace
            freq = float(self.get_count(go_id)) / float(self.get_total_count(namespace))
        except ZeroDivisionError:
            freq = 0

        return freq
