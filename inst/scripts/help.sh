#!/usr/bin/env bash
# PONG2 Help System
# KIR Genotype Imputation and Analysis Toolkit
# Version: 1.0.0 (update as needed)

set -euo pipefail

# ────────────────────────────────────────────────
# Color & Style Definitions
# ────────────────────────────────────────────────
BOLD=$(tput bold)
GREEN=$(tput setaf 2)
BLUE=$(tput setaf 4)
CYAN=$(tput setaf 6)
YELLOW=$(tput setaf 3)
RED=$(tput setaf 1)
WHITE=$(tput setaf 7)
RESET=$(tput sgr0)

# ────────────────────────────────────────────────
# Main Help Screen (overview of all commands)
# ────────────────────────────────────────────────
show_main_help() {
    cat << EOF
${BOLD}${GREEN}PONG2 ─ KIR Genotype Imputation & Training Toolkit${RESET}
${BLUE}Version:${RESET} 1.0.0
${BLUE}Description:${RESET} High-accuracy KIR allele imputation and model training from SNP data (chr19 locus)

${BOLD}${GREEN}USAGE${RESET}
  pong2 <command> [options]

${BOLD}${GREEN}AVAILABLE COMMANDS${RESET}
  ${CYAN}impute${RESET}     Predict KIR alleles from target genotype data
  ${CYAN}train${RESET}      Train a new KIR prediction model from reference data

${BOLD}${GREEN}GLOBAL OPTIONS${RESET}
  ${CYAN}-h, --help${RESET}     Show this help message (or command-specific help)
  ${CYAN}--version${RESET}      Show version information

${BOLD}${GREEN}GET COMMAND-SPECIFIC HELP${RESET}
  pong2 impute --help
  pong2 train --help

${BOLD}${GREEN}EXAMPLES${RESET}
  # Basic imputation
  pong2 impute -i data/target_chr19 -o results/impute -l KIR -a hg38

  # Imputation with missing SNP fill-in
  pong2 impute -i data/low_cov -o results/filled -l KIR -a hg38 --fill-missing --threads 16

  # Train a new model
  pong2 train -i data/ref_chr19 -k data/kir_calls.csv -o models/v2 -l KIR -a hg19

${BOLD}${YELLOW}IMPORTANT NOTES${RESET}
  • Input must be PLINK format (bed/bim/fam) covering the KIR region
  • Use --fill-missing for low SNP coverage (requires minimac4 4.1.6)
  • Always review the matching percentage report in the output directory

${BOLD}${BLUE}DOCUMENTATION & SUPPORT${RESET}
  GitHub:    https://github.com/NormanLabUCD/PONG2
  Issues:    https://github.com/NormanLabUCD/PONG2
  Citation:  Sadeeq et al. (2025). PONG 2.0: Allele imputation for KIR. Human Immunology.
EOF
}

# ────────────────────────────────────────────────
# Impute-specific detailed help
# ────────────────────────────────────────────────
show_impute_help() {
    cat << EOF
${BOLD}${GREEN}pong2 impute ─ KIR Allele Imputation${RESET}

${BOLD}${WHITE}DESCRIPTION${RESET}
  Predicts KIR genotypes from target SNP data using pre-trained models.
  Includes SNP quality checks and optional missing data imputation.

${BOLD}${WHITE}REQUIRED OPTIONS${RESET}
  ${CYAN}-i, --input FILE${RESET}     PLINK bed/bim/fam prefix (KIR region recommended)
  ${CYAN}-o, --output DIR${RESET}     Output directory (created if missing)
  ${CYAN}-l, --locus STR${RESET}      Target locus (currently only 'KIR')
  ${CYAN}-a, --assembly STR${RESET}   Genome build: 'hg19' or 'hg38'

${BOLD}${WHITE}OPTIONAL OPTIONS${RESET}
  ${CYAN}-t, --threads INT${RESET}    Number of CPU threads [default: 20]
  ${CYAN}--filter FLOAT${RESET}       Minimum SNP matching rate [default: 0.01]
  ${CYAN}-f, --force${RESET}          Proceed even if SNP matching rate is low
  ${CYAN}--fill-missing${RESET}       Impute missing SNPs using minimac4

${BOLD}${WHITE}EXAMPLES${RESET}
  Basic imputation:
    pong2 impute -i data/my_genotypes -o results/run1 -l KIR -a hg38

  With imputation fallback:
    pong2 impute -i low_coverage -o results/imputed -l KIR -a hg38 --fill-missing -t 32

  Force run:
    pong2 impute -i partial_data -o results/forced -l KIR -a hg19 -f

${BOLD}${WHITE}OUTPUT FILES${RESET}
  • predicted_kir_genotypes.csv    Main results: sample ID, predicted alleles, confidence
  • snp_matching_report.txt        SNP coverage & quality metrics
  • imputation_log.txt             Processing log
  • tmp/                           Temporary files (cleaned up by default)

${BOLD}${WHITE}TROUBLESHOOTING${RESET}
  • Low SNP matching → Try --fill-missing or pre-impute with Michigan Server
  • minimac4 not found → Install minimac4 4.1.6 and add to PATH
EOF
}

# ────────────────────────────────────────────────
# Train-specific detailed help
# ────────────────────────────────────────────────
show_train_help() {
    cat << EOF
${BOLD}${GREEN}pong2 train ─ Train a New KIR Prediction Model${RESET}

${BOLD}${WHITE}DESCRIPTION${RESET}
  Builds a new KIR prediction model from reference genotypes and known KIR calls.
  Useful for custom populations, updated reference panels, or locus-specific retraining.

${BOLD}${WHITE}REQUIRED OPTIONS${RESET}
  ${CYAN}-i, --bfile FILE${RESET}     PLINK bed/bim/fam prefix (reference genotypes)
  ${CYAN}-k, --kfile FILE${RESET}     CSV file with sample IDs and known KIR calls
  ${CYAN}-o, --output DIR${RESET}     Directory to save trained model files
  ${CYAN}-l, --locus STR${RESET}      Target locus (currently only 'KIR')
  ${CYAN}-a, --assembly STR${RESET}   Genome build: 'hg19' or 'hg38'

${BOLD}${WHITE}OPTIONAL OPTIONS${RESET}
  ${CYAN}-t, --threads INT${RESET}    Number of CPU threads [default: 20]
  ${CYAN}--filter FLOAT${RESET}       SNP quality filter threshold [default: 0.01]

${BOLD}${WHITE}EXAMPLES${RESET}
  Basic training:
    pong2 train -i data/reference_chr19 -k data/kir_calls.csv -o models/v2 -l KIR -a hg19

  With more threads and stricter filtering:
    pong2 train -i ref_panel -k calls.csv -o models/custom -l KIR -a hg38 -t 32 --filter 0.005

${BOLD}${WHITE}INPUT REQUIREMENTS${RESET}
  • --bfile: PLINK files covering the KIR region (chr19)
  • --kfile: CSV format, at minimum: sample_id,allele1,allele2
    (additional columns allowed and ignored)

${BOLD}${WHITE}OUTPUT FILES${RESET}
  • model.rds / model.RData           Trained prediction model
  • feature_importance.csv            Variable importance (if available)
  • training_log.txt                  Model training summary and metrics
  • tmp/                              Intermediate files (cleaned up by default)

${BOLD}${WHITE}TROUBLESHOOTING${RESET}
  • Missing KIR calls → Ensure --kfile has matching sample IDs
  • Low-quality SNPs → Increase --filter or pre-filter reference data
  • Training takes too long → Increase --threads or reduce dataset size
EOF
}

# ────────────────────────────────────────────────
# Dispatcher
# ────────────────────────────────────────────────
case "${1:-}" in
    "impute")
        show_main_help
        echo ""
        show_impute_help
        ;;
    "train")
        show_main_help
        echo ""
        show_train_help
        ;;
    "" | "-h" | "--help" | "help")
        show_main_help
        ;;
    "version")
        echo "PONG2 version 1.0.0"
        ;;
    *)
        echo -e "${RED}Error: Unknown command or option '$1'${RESET}" >&2
        echo -e "Use: pong2 --help" >&2
        exit 1
        ;;
esac

exit 0
