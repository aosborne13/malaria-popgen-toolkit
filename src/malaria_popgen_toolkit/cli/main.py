from malaria_popgen_toolkit.commands import missense_drugres_af

# ...
if args.command in ("missense-drugres-af", "run", "missense-af"):
    missense_drugres_af.run(
        vcf=args.vcf,
        ref_fasta=args.ref,
        gff3=args.gff3,
        metadata_path=args.metadata,
        outdir=args.outdir,
        min_dp=args.min_dp,
    )
