# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Independent R1 referee checks. No imports from the repository."""
import itertools as it, math, json, time
from collections import Counter, defaultdict
from fractions import Fraction
import numpy as np
import networkx as nx

rng = np.random.default_rng(20261103)
summary = {}

def norm(a): return float(np.linalg.svd(np.sqrt(a), compute_uv=False)[0])
def kernel(p, extreme=None):
    vertices=np.array([[1-2*p,p,p],[0,1-p,p],[0,p,1-p]])
    return vertices[extreme] if extreme is not None else rng.dirichlet([.15,.15,.15])@vertices

def make_chain(r,p,choices=None):
    # keys are separate row/column integers; kernels may depend on BOTH pasts.
    law={(0,0):1.0}; product=1.0
    for k in range(r):
        nxt={}; gamma=0
        for j,((x,y),w) in enumerate(sorted(law.items())):
            v=kernel(p, choices[k][j] if choices is not None else None)
            gamma=max(gamma,norm(np.array([[v[0],v[1]],[v[2],0.]])))
            for (a,b),z in zip([(0,0),(0,1),(1,0)],v):
                # Keep zero branches for complete extremal enumeration; they can only loosen gamma.
                nxt[(2*x+a,2*y+b)]=w*z
        law=nxt; product*=gamma
    a=np.zeros((2**r,2**r))
    for (x,y),w in law.items(): a[x,y]=w
    observed=norm(a); sharp=((1+math.sqrt(max(0,1-4*p*p)))/2)**(r/2)
    assert observed <= product+2e-12 and product <= sharp+2e-12
    return product-observed, observed, product

start=time.time(); worst=1.; tests=0
for p in [.001,.1,.25,.49,.5]:
    for labels in it.product(range(3),repeat=4):
        slack,_,_=make_chain(2,p,[[labels[0]],list(labels[1:])]); worst=min(worst,slack);tests+=1
for _ in range(2400):
    r=int(rng.integers(1,7)); p=float(rng.uniform(.0001,.5))
    slack,_,_=make_chain(r,p); worst=min(worst,slack);tests+=1
summary['chains']={'tests':tests,'min_product_minus_norm':worst,'seconds':time.time()-start}
print('CHAINS',json.dumps(summary['chains']),flush=True)

def edge_list(G): return sorted(tuple(sorted(e)) for e in G.edges())
def trees(G):
    E=edge_list(G); n=len(G); out=[]
    for ids in it.combinations(range(len(E)),n-1):
        par=list(range(n))
        def root(a):
            while par[a]!=a: a=par[a]
            return a
        good=True
        for i in ids:
            a,b=map(root,E[i])
            if a==b: good=False;break
            par[a]=b
        if good:out.append(sum(1<<i for i in ids))
    assert out
    return E,out

def cycles(G,lmax=4):
    E=edge_list(G); eid={e:i for i,e in enumerate(E)}; out={}
    for l in range(3,min(lmax,len(G))+1):
        for vs in it.combinations(range(len(G)),l):
            for tail in it.permutations(vs[1:]):
                seq=(vs[0],)+tail
                if seq[1]>seq[-1]:continue
                es=[tuple(sorted((seq[i],seq[(i+1)%l]))) for i in range(l)]
                if all(e in eid for e in es):out[sum(1<<eid[e] for e in es)]=frozenset(vs)
    return out

def stat(t,A,cs):
    B=~A; x=y=0
    for c in cs:
        x=2*x+int((t&(c&A))==(c&A)); y=2*y+int((t&(c&B))==(c&B))
    return x,y

def chain_from_counts(counts,r):
    a=np.zeros((2**r,2**r)); total=sum(counts.values())
    for xy,v in counts.items(): a[xy]=v/total
    product=1.; mincross=1.
    for k in range(r):
        blocks=defaultdict(lambda:np.zeros((2,2),dtype=np.int64))
        for (x,y),w in counts.items():
            blocks[(x>>(r-k),y>>(r-k))][((x>>(r-k-1))&1,(y>>(r-k-1))&1)]+=w
        gamma=0
        for block in blocks.values():
            b=block/block.sum(); assert b[1,1]==0
            mincross=min(mincross,b[0,1],b[1,0]); gamma=max(gamma,norm(b))
        product*=gamma
    return norm(a),product,mincross

atlas=[nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g() if 3<=len(g)<=5 and nx.is_connected(g)]
extras=[nx.circular_ladder_graph(3), nx.circular_ladder_graph(4), nx.cartesian_product(nx.complete_graph(4),nx.path_graph(2))]
extras=[nx.convert_node_labels_to_integers(g) for g in extras]
# A bridge between two disjoint triangles: earlier bits expose the outside but do not freeze the next cycle.
b=nx.disjoint_union(nx.complete_graph(3),nx.complete_graph(3)); b.add_edge(2,3);extras.append(b)
local_conditions=0; local_cycles=0; margin=Fraction(10**10); multi_cycle_tests=0;cut_tests=0;dp_tests=0
for gi,G in enumerate(atlas+extras):
    E,ts=trees(G); cs=cycles(G); m=len(E);D=max(dict(G.degree()).values())
    for c,vs in cs.items():
        local_cycles+=1
        incident=sum(1<<i for i,(a,b) in enumerate(E) if a in vs or b in vs)
        bins=defaultdict(list)
        for t in ts: bins[t&~incident].append(t)
        strongest_den=min(math.prod(G.degree(y) for y in vs if y!=root) for root in vs)
        for group in bins.values():
            pat=Counter(t&c for t in group)
            for i in range(m):
                if (c>>i)&1:
                    hits=pat[c^(1<<i)]; local_conditions+=1
                    assert hits*strongest_den>=len(group),(gi,c,group)
                    margin=min(margin,Fraction(hits*strongest_den,len(group)))
    if m<=10:
        cuts=[a for a in range(1,1<<m) if m/3<=a.bit_count()<=2*m/3]
    else:
        cuts=list({sum(1<<int(i) for i in rng.choice(m,int(rng.integers(math.ceil(m/3),math.floor(2*m/3)+1)),replace=False)) for _ in range(100)})
    for A in cuts:
        cut_tests+=1
        mixed=[c for c in cs if c&A and c&~A]
        packs=[]
        for c in mixed:
            packs.append([c])
            for d in mixed:
                if c<d and not cs[c]&cs[d]:packs.extend([[c,d],[d,c]])
        # Check every disjoint pair and its reverse on larger graphs, one single cycle otherwise.
        if any(len(p)==2 for p in packs): packs=[p for p in packs if len(p)==2]
        else:packs=packs[:1]
        if not packs:continue
        beta=None
        for pack in packs:
            counts=Counter(stat(t,A,pack) for t in ts)
            observed,product,mincross=chain_from_counts(counts,len(pack))
            p=D**(-(max(c.bit_count() for c in pack)-1))
            assert mincross+1e-14>=p
            assert observed <= product+2e-12 <= (1-p*p)**(len(pack)/2)+4e-12
            if len(pack)>1:multi_cycle_tests+=1
            # Avoid expensive UST SVD on the 16-edge graph except first ten cuts.
            if m<=10 or cut_tests%20==0:
                if beta is None:
                    xs={x:i for i,x in enumerate(sorted({t&A for t in ts}))};ys={y:i for i,y in enumerate(sorted({t&~A for t in ts}))}
                    mat=np.zeros((len(xs),len(ys)))
                    for t in ts:mat[xs[t&A],ys[t&~A]]=1/len(ts)
                    beta=norm(mat)
                assert beta <= observed+2e-12
                dp_tests+=1
summary['ust']={'graphs':len(atlas+extras),'cycles':local_cycles,'conditional_deletion_inequalities':local_conditions,'minimum_local_margin':str(margin),'cuts':cut_tests,'ordered_disjoint_pair_checks':multi_cycle_tests,'data_processing_checks':dp_tests}
print('UST',json.dumps(summary['ust']),flush=True)

def expansion(G):
    n=len(G);E=edge_list(G)
    return min(Fraction(sum(((s>>a)^(s>>b))&1 for a,b in E),s.bit_count()) for s in range(1,1<<n) if s.bit_count()<=n//2)

geometry_cuts=0; geometry_graphs=0
bases=[nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g() if 3<=len(g)<=4 and nx.is_connected(g)]+[nx.cycle_graph(5)]
for H in bases:
    hH=expansion(H); Delta=max(dict(H.degree()).values());n=len(H)
    for family in ['prism','square']:
        if family=='prism':
            G=nx.Graph();G.add_nodes_from(range(2*n));G.add_edges_from((u+k*n,v+k*n) for u,v in H.edges for k in (0,1));G.add_edges_from((u,u+n) for u in H)
            vcycles=[[(u,v),(v,v+n),(v+n,u+n),(u+n,u)] for u,v in H.edges]
            ell=4;kappa=Delta;h=min(Fraction(1),hH)/3
        else:
            G=nx.power(H,2);vcycles=[]
            for v in H:
                for u,w in it.combinations(H.neighbors(v),2):vcycles.append([(u,v),(v,w),(w,u)])
            ell=3;kappa=Fraction(3*Delta*(Delta-1),2);h=hH
        E=edge_list(G);ids={e:i for i,e in enumerate(E)};m=len(E);N=len(G);D=max(dict(G.degree()).values())
        cs={sum(1<<ids[tuple(sorted(e))] for e in c):frozenset(sum(([a,b] for a,b in c),[])) for c in vcycles}
        assert expansion(G)>=h
        for x in G:
            link=nx.Graph();link.add_nodes_from(i for i,(a,b) in enumerate(E) if x in (a,b))
            for c,vs in cs.items():
                if x in vs:
                    pair=[i for i,(a,b) in enumerate(E) if c>>i&1 and x in (a,b)];link.add_edge(*pair)
            assert nx.is_connected(link)
            assert sum(x in vs for vs in cs.values())<=kappa
        lower=min(h*N/(4*D*ell),Fraction(2*m,3)/(ell*D*(1+D/h)))
        for A in range(1,1<<m):
            if not m/3<=A.bit_count()<=2*m/3:continue
            mixed=[c for c in cs if c&A and c&~A];assert len(mixed)>=lower
            used=set();packed=0
            for c in mixed:
                if not used&cs[c]:packed+=1;used.update(cs[c])
            assert packed*ell*kappa>=len(mixed)
            geometry_cuts+=1
        geometry_graphs+=1
summary['geometry']={'graphs':geometry_graphs,'balanced_cuts':geometry_cuts}
# Uniform all spanning forests of a triangle: the unmodified local lemma is false.
summary['forest_counterexample']={'graph':'C3','forests':7,'fixed_deleted_edge_path_probability':'1/7','claimed_UST_local_lower_bound':'1/4'}
print('GEOMETRY',json.dumps(summary['geometry']),flush=True)
print('SUMMARY',json.dumps(summary,sort_keys=True),flush=True)
