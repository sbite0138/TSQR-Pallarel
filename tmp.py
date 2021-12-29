
import random


def solve(block_num, proc_num, have1, have2, need1, need2):
    have1, have2, need1, need2 = need1, need2, have1, have2

    send1 = -1
    send2 = -1
    recv1 = -1
    recv2 = -1

    recv1 = need1 % proc_num
    recv2 = -1 if need2 == -1 else need2 % proc_num

    t = proc_num - (block_num % proc_num)
    send1 = have1
    if (send1 >= t):
        send1 = (have1-t)//2+t

    if (have2 != -1):
        send2 = have2
        if (send2 >= t):
            send2 = (have2-t)//2+t
    send1, send2, recv1, recv2 = recv1, recv2, send1, send2
    return send1, send2, recv1, recv2


def calc_have_need(rank, block_num, proc_num):
    have1 = -1
    have2 = -1
    need1 = -1
    need2 = -1
    if (rank < block_num % proc_num):
        have1 = rank
        have2 = proc_num + rank
    else:
        have1 = rank
    t = proc_num - (block_num % proc_num)
    if (rank < t):
        need1 = rank
    else:
        need1 = 2*(rank - t) + t
        need2 = need1+1
    have1, have2, need1, need2 = need1, need2, have1, have2
    return have1, have2, need1, need2


def judge(block_num, proc_num):
    have_table = [[] for _ in range(proc_num)]
    need_table = [[] for _ in range(proc_num)]
    need_all = []
    have_all = []

    for rank in range(proc_num):
        have1, have2, need1, need2 = calc_have_need(rank, block_num, proc_num)
        assert(have1 != -1)
        assert(have1 not in have_all)
        have_all.append(have1)
        if (have2 != -1):
            assert(have2 not in have_all)
            have_all.append(have2)

        assert(need1 != -1)
        assert(need1 not in need_all)
        need_all.append(need1)
        if (need2 != -1):
            assert(need2 not in need_all)
            need_all.append(need2)
        have_table[rank].append(have1)
        have_table[rank].append(have2)
        need_table[rank].append(need1)
        need_table[rank].append(need2)
    print(have_all, block_num)
    print(need_all, block_num)
    print(len(have_all), block_num)
    assert(len(have_all) == block_num)
    assert(len(need_all) == block_num)
    for i in range(block_num):
        assert(i in have_all)
        assert(i in need_all)

    for rank in range(proc_num):
        have1, have2, need1, need2 = calc_have_need(rank, block_num, proc_num)
        send1, send2, recv1, recv2 = solve(
            block_num, proc_num, have1, have2, need1, need2)
        assert(send1 == -1 or have1 in need_table[send1])
        assert(send2 == -1 or have2 in need_table[send2])
        assert(recv1 == -1 or need1 in have_table[recv1])
        assert(recv2 == -1 or need2 in have_table[recv2])


for i in range(100):
    block_num = random.randint(1, 1000)
    proc_num = random.randint(1, 1000)
    while (not(block_num >= proc_num and 2*proc_num > block_num)):
        block_num = random.randint(1, 1000)
        proc_num = random.randint(1, 1000)
    print(i, block_num, proc_num)
    judge(block_num, proc_num)
