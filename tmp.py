


m=10
proc_num=4
m_part = m // proc_num
block_num = (m + m_part - 1) // m_part
block_num=7
for rank in range(proc_num):
    have1 = -1
    have2 = -1
    need1 = -1
    need2 = -1
    if (rank < block_num % proc_num):
        need1 = 2 * rank
        need2 = 2 * rank + 1
        have1 = rank
        have2 = proc_num + rank
    else:
        need1 = rank + (block_num % proc_num)
        have1 = rank
    print(f"rank={rank}, have1={have1}, have2={have2}, need1={need1}, need2={need2}")