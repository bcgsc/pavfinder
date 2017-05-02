from intspan import intspan
from sets import Set
from alignment import reverse_complement, compare_chr
from adjacency import Adjacency

def find_chimera(aligns, query_seq=None, max_splits=3, debug=False):
    path = find_paths(aligns, min_coverage=0.01, from_edge=0.8, debug=debug)
    adjs = []
    if path:
	chimeric_aligns = [aligns[idx] for idx in path]
	for i in range(len(chimeric_aligns) - 1):
	    adj = call_event(chimeric_aligns[i], chimeric_aligns[i + 1], query_seq=query_seq, no_sort=True)
	    if adj is not None:
		if debug:
		    print 'chimera1', chimeric_aligns[i].query, chimeric_aligns[i].target, \
		          chimeric_aligns[i].tstart, chimeric_aligns[i].tend, \
		          chimeric_aligns[i].qstart, chimeric_aligns[i].qend, \
		          chimeric_aligns[i].strand
		    print 'chimera2', chimeric_aligns[i + 1].query, chimeric_aligns[i + 1].target, \
		          chimeric_aligns[i + 1].tstart, chimeric_aligns[i + 1].tend, \
		          chimeric_aligns[i + 1].qstart, chimeric_aligns[i + 1].qend, \
		          chimeric_aligns[i + 1].strand
		    
		# check for inverted duplication
		check_inv_dup(adj, (chimeric_aligns[i], chimeric_aligns[i + 1]))

		adjs.append(adj)

    return adjs

def check_inv_dup(adj, aligns):
    if adj.rearrangement == 'inv':
	target_span_before_bp = intspan('%s-%s' % (aligns[0].tstart, aligns[0].tend))
	target_span_after_bp = intspan('%s-%s' % (aligns[1].tstart, aligns[1].tend))
	if target_span_after_bp < target_span_before_bp:
	    adj.rearrangement = 'inv-dup'
	    
	    # reverse breakpoint and orientation to make it same as a dup
	    
	    if adj.target_breaks[1] == aligns[1].tstart:
		adj.target_breaks[1] = aligns[1].tend
	    else:
		adj.target_breaks[1] = aligns[1].tstart
		
	    if adj.orients[1] == 'L':
		adj.orients[1] = 'R'
	    else:
		adj.orients[1] = 'L'
		
def call_event(align1, align2, query_seq=None, no_sort=False, max_inv_target_olap=30000, debug=False):
    """Curates adj based on info given by primary_aligns alignments
    
    Args:
        align1: First Alignment object
	align2: Second Alignment object
	homol_seq: (str) Microhomology sequence in contig
	homol_coords: (tuple) Start and end coordinates of microhomology sequence in contig
	novel_seq: (str) Untemplated sequence at breakpoint
	query_seq: (str) Query sequence
	
    Returns:
	Adjacency object
    """    
    # figure out breakpoints using query positions
    target_breaks = [None, None]
    orients = [None, None]
    query_breaks = [None, None]
    homol_seq = None
    homol_seq_coords = None
    novel_seq = None
    novel_seq_coords = None

    align1_tpos = (align1.tstart, align1.tend) if align1.strand == '+' else (align1.tend, align1.tstart)
    align2_tpos = (align2.tstart, align2.tend) if align2.strand == '+' else (align2.tend, align2.tstart)
    if align1.qstart < align2.qstart:
	aligns = [align1, align2]
	target_breaks[0] = align1.tend if align1.strand == '+' else align1.tstart
	orients[0] = 'L' if max(align1.tstart, align1.tend) == target_breaks[0] else 'R'
	target_breaks[1] = align2.tstart if align2.strand == '+' else align2.tend
	orients[1] = 'L' if max(align2.tstart, align2.tend) == target_breaks[1] else 'R'
	query_breaks = [align1.qend, align2.qstart]
	
    else:
	aligns = [align2, align1]
	target_breaks[0] = align2.tend if align2.strand == '+' else align2.tstart
	orients[0] = 'L' if max(align2.tstart, align2.tend) == target_breaks[0] else 'R'
	target_breaks[1] = align1.tstart if align1.strand == '+' else align1.tend
	orients[1] = 'L' if max(align1.tstart, align1.tend) == target_breaks[1] else 'R'
	query_breaks = [align2.qend, align1.qstart]
	
    if not no_sort:
	if (aligns[0].target != aligns[1].target and compare_chr(aligns[0].target, aligns[1].target) > 0) or\
	   (aligns[0].target == aligns[1].target and target_breaks[0] > target_breaks[1]):
	    aligns.reverse()
	    target_breaks.reverse()
	    orients.reverse()
	    	
    rearrangement = None
    if aligns[0].target != aligns[1].target:
	rearrangement = 'trl'
    elif orients[0] == orients[1]:
	span1 = intspan('%s-%s' % (aligns[0].tstart, aligns[0].tend))
	span2 = intspan('%s-%s' % (aligns[1].tstart, aligns[1].tend))
	olap = span1 & span2
	if len(olap) <= max_inv_target_olap:
	    rearrangement = 'inv'
	else:
	    print '%s:potential inv disallowed - target overlap %d bigger than %s' % (aligns[0].query,
	                                                                           len(olap),
	                                                                           max_inv_target_olap)
    elif orients[0] == 'L' and orients[1] == 'L':
	rearrangement = 'inv'
    elif orients[0] == 'L' and orients[1] == 'R':
	if target_breaks[0] < target_breaks[1]:
	    if target_breaks[0] + 1 == target_breaks[1]:
		# deletion of tandem duplicaton
		if query_breaks[0] >= query_breaks[1]:
		    rearrangement = 'del'
		    target_breaks = [target_breaks[1] + 1, target_breaks[0] + (query_breaks[0] - query_breaks[1] + 1)]
		else:
		    rearrangement = 'ins'
	    else:
		# deletion with or without microhology
		rearrangement = 'del'
			    
	elif target_breaks[0] > target_breaks[1]:
	    rearrangement = 'dup'

	else:
	    if query_breaks[0] < query_breaks[1]:
		rearrangement = 'ins'
	    else:
		# deletion of tandem duplicaton
		rearrangement = 'del'
		target_breaks = [target_breaks[1] + 1, target_breaks[0] + (query_breaks[0] - query_breaks[1] + 1)]
	
    elif orients[0] == 'R' and orients[1] == 'R':
	rearrangement = 'inv'
    elif orients[0] == 'R' and orients[1] == 'L':
	if target_breaks[0] == target_breaks[1]:
	    rearrangement = 'ins'
	elif target_breaks[0] < target_breaks[1]:
	    rearrangement = 'dup'
	else:
	    rearrangement = 'del'
	    
    # novel seq
    if query_seq is not None and query_breaks[1] - query_breaks[0] > 1:
	novel_seq = query_seq[query_breaks[0] : query_breaks[1] - 1]
	if aligns[0].strand == '-':
	    novel_seq = reverse_complement(novel_seq)
	    novel_seq_coords = (query_breaks[0] + 1, query_breaks[1] - 1)

    # homol seq
    if query_seq is not None and query_breaks[0] >= query_breaks[1]:
	homol_seq_coords = [query_breaks[1], query_breaks[0]]
	homol_seq = query_seq[query_breaks[1] - 1 : query_breaks[0]]
	if aligns[0].strand == '-':
	    homol_seq = reverse_complement(homol_seq)
	    homol_seq_coords = (query_breaks[1], query_breaks[0])
			
    adj = None
    if rearrangement is not None:
	adj = Adjacency(align1.query,
	                (aligns[0].target, aligns[1].target),
	                query_breaks,
	                target_breaks,
	                rearrangement = rearrangement,
	                orients = orients,
	                homol_seq = homol_seq,
	                homol_seq_coords = homol_seq_coords,
	                novel_seq = novel_seq,
	                novel_seq_coords = novel_seq_coords,
	                )
	    
    elif debug:
	sys.stdout.write("cannot figure out event of primary_aligns alignment contig:%s targets:%s,%s orients:%s breaks:%s query_breaks:%s\n" % (aligns[0].query,
	                                                                                                                                         aligns[0].target,
	                                                                                                                                         aligns[1].target,
	                                                                                                                                         orients,
	                                                                                                                                         breaks,
	                                                                                                                                         query_breaks))
    return adj

    
def find_paths(aligns, min_coverage=None, use_end_to_end=True, get_all=False, max_nodes=500, max_paths=5, max_ends=50, no_trim=[], same_target=None, from_edge=0.3, debug=False):
    def _find_end_points():
	starts = []
	ends = []
	from_start = max(1, int(from_edge * aligns[0].query_len))
	from_end = aligns[0].query_len - from_start + 1
	for i in range(len(aligns)):
	    if aligns[i].qstart <= from_start:
		try:
		    starts.append(i)
		except:
		    starts = [i]
	    if aligns[i].qend >= from_end:
		try:
		    ends.append(i)
		except:
		    ends = [i]
		
	return starts, ends
    
    def _trim_aligns(max_mappings=100):
	regions = {}
	for i in range(len(aligns)):
	    if i in no_trim:
		continue
	    
	    key = '%s-%s:%s' % (aligns[i].qstart, aligns[i].qend, aligns[i].strand)
	    try:
		regions[key].append(i)
	    except:
		regions[key] = [i]
		
	extra = []
	for key in regions.keys():
	    if len(regions[key]) > max_mappings:
		for i in regions[key][max_mappings:]:
		    extra.append(i)
		
	if extra:
	    new_aligns = []
	    for i in range(len(aligns)):
		if not i in extra:
		    new_aligns.append(aligns[i])
		    
	    return new_aligns
	
	return None

    def _construct_graph(max_links=100, max_olap=0.2, max_gap=0.2):
	"""
	Alignments must be like to be an edge
	   ----> j
	----> i
	"""
	graph = {}
	for i in range(len(aligns)):
	    len_i = aligns[i].qend - aligns[i].qstart + 1
	    for j in range(len(aligns)):
		if i != j:
		    len_j = aligns[j].qend - aligns[j].qstart + 1
		    olap = max(0, aligns[i].qend - aligns[j].qstart + 1)
		    gap = max(0, aligns[j].qstart - aligns[i].qend - 1)
		    # only allows gap = 0, gap > 0, or partial overlap
		    if aligns[j].qend > aligns[i].qend and aligns[j].qstart > aligns[i].qstart:
			if float(olap)/len_i >= max_olap or\
			   float(olap)/len_j >= max_olap or\
			   float(gap)/len_i >= max_gap or\
			   float(gap)/len_j >= max_gap:
			    continue
			try:
			    graph[i].append(j)
			except:
			    graph[i] = [j]
			    
	return graph

    
    def _bfs(graph, start, end):
	"""Breath-first search
	http://stackoverflow.com/questions/8922060/breadth-first-search-trace-path
	"""
	paths = []
	# maintain a queue of paths
	queue = []
	# push the first path into the queue
	queue.append([start])
	while queue:
	    # get the first path from the queue
	    path = queue.pop(0)
	    # get the last node from the path
	    node = path[-1]
	    # path found
	    if node == end:
		# captures all paths
		paths.append(path)
		# this will return shortest path
		#return path
	    # enumerate all adjacent nodes, construct a new path and push it into the queue
	    for adjacent in graph.get(node, []):
		new_path = list(path)
		new_path.append(adjacent)
		queue.append(new_path)
		
	return paths
    
    def _coverage(path):
	spans = [intspan('%d-%d' % (aligns[i].qstart, aligns[i].qend)) for i in path]
	covered = spans[0]
	overlaps = []
	for i in range(1, len(spans)):
	    covered = covered.union(spans[i])
	    
	    overlap = spans[i - 1].intersection(spans[i])
	    if len(overlap) > 0:
		overlaps.append(overlap)
	    
	return covered, overlaps
    
    def _screen(covered, overlaps):
	passed = True
	if min_coverage is not None:
	    if not use_end_to_end:
		coverage = len(covered) / float(aligns[0].query_len)
	    else:
		coverage = (max(covered) - min(covered) + 1) / float(aligns[0].query_len)

	    if coverage < min_coverage:
		passed = False
		
	return passed
    
    def _pick_best(path_info):
	best = {'index': None, 'covered': None, 'overlapped': None}
	for i in sorted(path_info.keys()):
	    covered = len(path_info[i]['covered'])
	    if not path_info[i]['overlaps']:
		overlapped = 0
	    else:
		overlapped = sum([len(olap) for olap in path_info[i]['overlaps']])

	    if best['index'] is None or\
	       covered > best['covered'] or\
	       overlapped < best['overlapped']:
		best['index'] = i
		best['covered'] = covered
		best['overlapped'] = overlapped
	
	return best['index']
    
	    
    # check promiscuity
    trimmed_aligns = _trim_aligns()
    if trimmed_aligns:
	aligns = trimmed_aligns
	
    # same target
    if same_target:
	aligns = [align for align in aligns if align.target == same_target]
    
    if not aligns:
	return []
    
    # construct graph
    graph = _construct_graph()
	    
    if len(graph.keys()) > max_nodes:
	if debug:
	    sys.stdout.write('%s: too many nodes(%d) to construct path(max:%d)\n' % (aligns[0].query, len(graph.keys()), max_nodes))
	return []
    
    # get end points
    starts, ends = _find_end_points()
    if not starts or not ends:
	return []
    
    if len(starts) > max_ends:
	starts = starts[:max_ends]
    if len(ends) > max_ends:
	ends = ends[:max_ends]
    
    # find paths
    paths = []
    for start in starts:
	if max_paths and len(paths) >= max_paths:
	    break
	for end in ends:
	    paths.extend(_bfs(graph, start, end))
	    if max_paths and len(paths) >= max_paths:
		break
	    	
    # screen
    path_info = {}
    for i in range(len(paths)):
	if len(paths[i]) == 1:
	    continue
	covered, overlaps = _coverage(paths[i])
	if _screen(covered, overlaps):
	    chroms = Set([aligns[j].target for j in paths[i]])
	    path_info[i] = {'path': paths[i], 'covered': covered, 'overlaps': overlaps, 'chroms':chroms}
	    	
    if path_info:
	if get_all:
	    return [paths[i] for i in path_info.keys()]
	else:
	    return paths[_pick_best(path_info)]
    else:
	return []
