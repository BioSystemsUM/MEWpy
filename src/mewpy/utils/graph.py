import networkx as nx
from mewpy.simulation.simulation import Simulator
from mewpy.simulation import get_simulator
import numpy as np
import math


def create_metabolic_graph(model,directed=True,reactions=None,remove=[],edges_labels=False):
    """ Creates a metabolic graph

    
    :param model: A model or a model containter
    :param (bool) directed: Defines if the graph to be directed or undirected. Defaults to True.
    :param (list) reactions: List of reactions to be included in the graph. Defaults to None, in which\
        all reactions are included.
    :param list remove: list os metabolites not to be included. May be used to remove cofactores such as ATP/ADP, NAD(P)(H), and acetyl-CoA/CoA.    
    :param (bool) edges_labels: Adds a reversabily label to edges. Defaults to False.
    :returns: A networkx graph of the metabolic network.
    """

    if not isinstance(model,Simulator):
        container = get_simulator(model)
    else:
        container = model

    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    if not reactions:
        reactions = container.reactions

    reactions = list(set(reactions)-set(remove))    

    for r in reactions:
            G.add_node(r, label=r, node_class="reaction", node_id=r)
    
    for r in reactions:
        the_metabolites = container.get_reaction_metabolites(r)
        for m in the_metabolites:
            if m in remove:
                continue
            if m not in G.nodes:
                G.add_node(m, label=m, node_class="metabolite", node_id=m)
            # evaluating if the metabolite has been defined as a reactant or product
            if the_metabolites[m]<0:
                (tail,head) = (m,r)
            elif the_metabolites[m]>0:
                (tail,head) = (r,m)

            # adding an arc between a metabolite and a reactions
            G.add_edge(tail,head)
            label = 'irrev'
            lb,_ = container.get_reaction_bounds(r)

            if lb<0:
                G.add_edge(head,tail)
                label = 'rev'

            if edges_labels:
                G[tail][head]['label'] = label

            G[tail][head]['reversible'] = lb<0

    return G


def shortest_distance(model,reaction,reactions=None, remove=[]):
    """ Returns the unweighted shortest path distance from a list of reactions to a reaction.
    Distances are the number of required reactions. If there is no pathway between the reactions the distance is inf·

    :param model: A model or a Simulator instance.
    :param str reaction: target reaction.
    :param list reactions: List os source reactions. Defaults to None, in which case all model reactions are considered.
    :param list remove: List os metabolites not to be included. May be used to remove path that include \
        cofactores such as ATP/ADP, NAD(P)(H), and acetyl-CoA/CoA.    
    :returns: A dictionary of distances.
    
    """
    if not isinstance(model,Simulator):
        container = get_simulator(model)
    else:
        container = model

    rxns = reactions if reactions else container.reactions
    if reaction not in rxns:
        rxns.append(reaction)

    G = create_metabolic_graph(container,reactions=rxns,remove=remove)    
    sp = dict(nx.single_target_shortest_path_length(G,reaction))
    

    distances = {}
    for rxn in rxns:
        if rxn in sp:
            distances[rxn]=sp[rxn]//2
        else:
            distances[rxn]= np.inf
    return distances


def probabilistic_reaction_targets(model,product,targets,factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, reactions closer to the product are repeated more often in the new target list.
    Moreover, reactions from which there is no path (pathway or cofactors usage) to the product are removed.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target reactions.
    :param int factor: Maximum number of repetitions, also the distance after which all reactions are\
        considered with equal probability. Defaults to 10.
    :returns: A probabilistic target list.
    
    """
    distances = shortest_distance(model,product,targets)
    prob_targets=[]
    for t in targets:
        if distances[t]==np.inf or distances[t]==0:
            continue
        else:
            coef = math.ceil(1/distances[t]*factor)
        x = [t]*coef
        prob_targets.extend(x)        
    return prob_targets


def probabilistic_gene_targets(model,product,targets,factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, genes on GPRs of reactions closer to the product are repeated more
    often in the new target list.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target genes.
    :param int factor: Maximum number of repetitions. Defaults to 10.
    :returns: A probabilistic target list.
    """
    if not isinstance(model,Simulator):
        container = get_simulator(model)
    else:
        container = model

    # Reaction targets
    if not targets:
        genes = container.genes
    else:
        genes = targets

    rxns = container.get_reactions_for_genes(genes)
    rxn_distances = shortest_distance(model,product,rxns)

    # genes distances are the maximum of all reaction 
    # distances that they catalyse.

    prob_targets =[]

    for gene in genes:
        rxs = container.get_reactions_for_genes([gene])
        dd = [d for r,d in rxn_distances.items() if r in rxs]
        d =max(dd)
        if d==np.inf:
            coef = 1        
        elif d==0:
            continue
        else:
            coef = math.ceil(1/d*factor)
        x = [gene]*coef
        prob_targets.extend(x)        

    return prob_targets

    
def probabilistic_protein_targets(model,product,targets,factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, proteins used in reactions closer to the product are repeated 
    more often in the new target list.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target genes.
    :param int factor: Maximum number of repetitions. Defaults to 10.
    :returns: A probabilistic target list.
    """
    if not isinstance(model,Simulator):
        container = get_simulator(model)
    else:
        container = model
    raise NotImplementedError