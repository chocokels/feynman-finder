# October 9, 2020 - Kelsey Bates

import networkx as nx

def FeynmanFinderFunc(ket,bra,graph,eta,heterodyne=False):
    '''Recursive function to find Feynman diagrams of a system.
    
    Args:
        ket:        Name of starting ket
        bra:        Name of starting bra
        graph:      Directed graph of energy levels (using networkx). Directions
                    should point from ground to excited states. The names of the
                    nodes should include the names of the starting bra and ket.
        eta:        List of signs of signal pathways. E.g. [-1,1,1,-1] for S1.
        heterodyne: True for photoluminescence, False for heterodyne. Note that
                    heterodyne requires last eta to be -1 to work properly.
    
    Returns:
        feyn: A list of dicts. Each item in the list corresponds to a single
              Feynman diagram. The dict has four keys.
                  ket:  List of the kets in the Feynman diagram (left side
                        states).
                  bra:  List of the bras in the Feynman diagram (right side
                        states).
                  sign: A list of the signs of the signal pathway (essentially
                        just a copy of eta).
                  side: A list of the side of the interaction on the Feynman
                        diagram. 1 is left and -1 is right.
    '''
    feyn = []
    
    if not eta:
        if ket == bra and (heterodyne or [n for n in graph.predecessors(ket)]):
            feyn = [dict(ket = [ket],
                         bra = [bra],
                         sign = [],
                         side = [])]
    else:
        sign = eta[0]
        side = 1
        if sign == 1:
            states = [n for n in graph.successors(ket)]
        else:
            states = [n for n in graph.predecessors(ket)]
        for state in states:
            news = FeynmanFinderFunc(state,bra,graph,eta[1:],heterodyne);
            for new in news:
                feyn.append(dict(ket = [ket] + new['ket'],
                                 bra = [bra] + new['bra'],
                                 sign = [sign] + new['sign'],
                                 side = [side] + new['side']))
        side = -1
        if sign == -1:
            states = [n for n in graph.successors(bra)]
        else:
            states = [n for n in graph.predecessors(bra)]
        if not heterodyne or (len(eta) > 1):
            for state in states:
                news = FeynmanFinderFunc(ket,state,graph,eta[1:],heterodyne);
                for new in news:
                    feyn.append(dict(ket = [ket] + new['ket'],
                                     bra = [bra] + new['bra'],
                                     sign = [sign] + new['sign'],
                                     side = [side] + new['side']))
    return feyn

def PrintFeynman(feyn):
    '''Prints Feynman diagrams in a human-readable form.

    Not identical to a hand-drawn diagram, as vectors pointing away from the
    diagram will be a row off. I optimized simplicity over perfection.

    Args:
        Output of FeynmanFinderFunc, i.e. a list of dicts.
    '''
    for curr in feyn:
        print(' |' + curr['ket'][-1] + '><' + curr['bra'][-1] + '| ')
        for i in [-n -1 for n in range(len(feyn[1]['sign']))]:
            if curr['sign'][i] == 1:
                symb = '/'
            elif curr['sign'][i] == -1:
                symb = '\\'
            left = ' '
            right = ' '
            if curr['side'][i] == 1:
                left = symb
            elif curr['side'][i] == -1:
                right = symb
            print(left + '|' + curr['ket'][i - 1] + '><' + curr['bra'][i - 1] + '|' + right)
        print(' ')
