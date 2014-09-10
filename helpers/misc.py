'''
extra functions
'''
import sys

def get_state(p, q, op, oq):
    ''' computes state'''

    if op == 0 and oq == 0:
        return 0
    if op == 0 and oq == 1:
        return 1
    if op == 1 and oq == 0:
        return 2
    if op == 1 and oq == 1:
        return 3


def get_dist(p, q, state, width1, width2, left1, right1, left2, right2, ins_size):
    ''' given a state, determines distance between reads implied by that state  '''

    if state == 0:
        dr1 = width1 - left1
        dr2 = right2
    elif state == 1:
        dr1 = width1 - left1
        dr2 = width2 - left2
    elif state == 2:
        dr1 = right1
        dr2 = right2
    else:
        dr1 = right1
        dr2 = width2 - left2

    # finalize dist.
    dist = ins_size - dr1 - dr2

    # return dist.
    return dist

def get_orien(sam1, sam2, pair_mode):
    ''' gets orientation from SAM object for pairs'''


    # turn into int.
    if sam1.OFLAG == "0":
        o1 = 0
    else:
        o1 = 1

    if sam2.OFLAG == "0":
        o2 = 0
    else:
        o2 = 1

    # flip to ff pairing.
    if pair_mode == 1:
        o2 = 1 - o2
    if pair_mode == 2:
        o1 = 1 - o1

    # return them.
    return o1, o2

