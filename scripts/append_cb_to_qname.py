"""
append_cb_to_qname.py
将 BAM 文件中每条 read 的 CB tag 追加到 read name 末尾
输入:  input.bam
输出:  output.bam  (read name 格式: original_name|CB序列)

用法: python append_cb_to_qname.py input.bam output.bam
"""

import pysam
import sys

def append_cb_to_qname(in_bam, out_bam):
    n_total = 0
    n_with_cb = 0
    n_no_cb = 0

    with pysam.AlignmentFile(in_bam, 'rb') as fin, \
         pysam.AlignmentFile(out_bam, 'wb', header=fin.header) as fout:

        for read in fin:
            n_total += 1
            try:
                cb = read.get_tag('CB')
                read.query_name = read.query_name + '|' + cb
                n_with_cb += 1
            except KeyError:
                # 没有 CB tag，read name 不变
                n_no_cb += 1
            fout.write(read)

            if n_total % 5000000 == 0:
                sys.stderr.write(f'[INFO] Processed {n_total:,} reads, '
                                 f'CB found: {n_with_cb:,}, '
                                 f'no CB: {n_no_cb:,}\n')

    sys.stderr.write(f'[INFO] Done. Total: {n_total:,}, '
                     f'CB appended: {n_with_cb:,}, '
                     f'no CB (unchanged): {n_no_cb:,}\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: python append_cb_to_qname.py input.bam output.bam\n')
        sys.exit(1)
    append_cb_to_qname(sys.argv[1], sys.argv[2])
