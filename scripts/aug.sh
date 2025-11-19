#!/bin/sh

module load conda
mamba activate liftnew

augustus \
	--species=subCanephora \
	--UTR=on \
	--extrinsicCfgFile=/quobyte/bcmeyersgrp/tell/testingCofSM/subCanephora/training/subCanephora.extrinsic.M.RM.E.W.P.cfg \
	--allow_hinted_splicesites=atac \
	/quobyte/bcmeyersgrp/tell/testingCofSM/subCanephora/subCanephora_Cara_hifi_final.fa \
	--codingseq=on \
	--protein=on \
	--outfile=/quobyte/bcmeyersgrp/tell/testingCofSM/subCanephora/training/predicted.UTR.C.subCanephora.gff \
	--progress=true \
	--genemodel=complete \
	--hintsfile=/quobyte/bcmeyersgrp/tell/testingCofSM/subCanephora/subCanephora_lifted_hints.gff3 \
	--errfile=/quobyte/bcmeyersgrp/tell/testingCofSM/subCanephora/training/errors_subCanephora.UTR.err
