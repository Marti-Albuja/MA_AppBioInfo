SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

##--------------------------------------------##

usage:
	@echo "Usage: make <target>"
	@echo "Targets:"
	@echo "  simulate: Simulate RNA-seq reads"
	@echo "  identify: Identify differentially expressed genes"
	@echo "  pca: Perform PCA analysis"
	@echo "  heatmap: Generate heatmap of differentially expressed genes"
	@echo "  all: Run all targets"
	@echo "  clean: Remove all generated files"

simulate:
	# Simulate RNA-seq reads
	src/r/simulate_counts.r -r 5 -n 30000

identify:
	# Identify differentially expressed genes
	Rscript src/r/edger.r

	# Compare resulting and expected counts of differentially expressed genes 
	Rscript src/r/evaluate_results.r -a counts.csv -b edger.csv

pca:
	# Perform PCA analysis
	src/r/plot_pca.r -c edger.csv

heatmap:
	# Generate heatmap of differentially expressed genes
	src/r/plot_heatmap.r -c edger.csv

all: simulate identify pca heatmap

clean:
	rm -rf *.csv *.png *.pdf 

.PHONY: usage simulate identify pca heatmap all clean





