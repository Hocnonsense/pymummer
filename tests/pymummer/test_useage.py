# -*- coding: utf-8 -*-
"""
* @Date: 2024-08-12 17:25:29
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2024-08-15 20:09:04
* @FilePath: /pymummer/tests/pymummer/test_useage.py
* @Description:
"""
# """

from io import StringIO

from pymummer import delta, usage, flatten
from tests import Path, temp_output, test_files, test_temp


def test_doc():
    import doctest

    doctest.testmod(usage, raise_on_error=True)


delta_file = test_files / "MarsFilter2.delta"


def test_report_indel_looong():
    with StringIO() as buf:
        assert (1, 0) == usage.report_indel_looong(
            delta.Delta(delta_file, {}), stdout=buf
        )
        buf.seek(0)
        assert buf.readlines() == [
            "DeltaContig2(NODE_1564_length_766_cov_111365.326301[0:238](+), NZ_CP008888.1[3367:3605](-) ..1)\n",
            "[ masked ]TGACCCTGAACCTCAGGACTCTGAGCCTCAGGA NODE_1564_length_766_cov_111365.326301 [0:238](+)\n",
            "   205    |||||AG||G|||||AA|||||||||||A|||| 0 : 0\n",
            "[bp same ]TGACCAGGAGCCTCAAAACTCTGAGCCTAAGGA NZ_CP008888.1 [3367:3605](-)\n",
            "CCTGAGTCTGAC[ masked ]AGGTT-AGGCT[ masked ] NODE_1564_length_766_cov_111365.326301 [197:766](+)\n",
            "||||||C|||||    46    |||||+|||||   501     1 : 0\n",
            "CCTGAGCCTGAC[bp same ]AGGTTGAGGCT[bp same ] NZ_CP008888.1 [2970:3540](-)\n",
            "try to merge above alignments:\n",
            "[ masked ]TGACCCTGAACCTCAGGACTCTGAGCCTCAGGAGACTACGGTCGAGGGGAAGGTTAGGCTAAGGACCTTGAAACCTCAAAGTCTGAGCCTGTCCAAATCGAGCAGAAGGTTGAGGCGTCTAAACCTGAGTCTGACCAGGAGCCTCAAAACTCTGAGCCTAAGGATTCTCCCGATGATGTTCCTGAGTCTGACCGACTC[ masked ] NODE_1564_length_766_cov_111365.326301 [0:766](+)\n",
            "   205    |||||---|------|||----|--|||||--|-|--||--||----------|--|---|-|--|||--||--------|---||---|-||----||---|------------||---|----|||-|--|----------|------||----|-|--|-|-------|-|------|-|||--|---|-|||||   363     0 : -131\n",
            "[bp same ]TGACC---A------GGA----G--CCTCA--A-A--AC--TC----------T--G---A-G--CCT--AA--------G---GA---T-TC----TC---C------------CG---A----TGA-T--G----------T------TC----C-T--G-A-------G-T------C-TGA--C---C-GACTC[bp same ] NZ_CP008888.1 [2970:3605](-)\n",
        ]
