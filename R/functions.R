library(DESeq2)
library(openxlsx)

pseudobulk_DEseq2 = function(dds, contrast.in = c("condition", "SA", "PBS")) {
    
    dds <- DESeq(dds)
    
    contrast = contrast.in
    
    res <- results(dds, 
                   contrast = contrast,
                   alpha = 0.05)
    
    
    coef.ind = which(resultsNames(dds) == 
                         paste0( c(contrast.in[1], contrast.in[2], "vs", contrast.in[3]),
                                 collapse = "_") )
    
    res <- lfcShrink(dds = dds, 
                     # contrast =  contrast,
                     coef = coef.ind,
                     res=res)
    
    res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        arrange(pvalue, abs(log2FoldChange)) %>% 
        as_tibble()
    
    return(res_tbl)
}


plot_volcano = function(res_tbl, title, padj.thresh = 0.05, lfc.thresh = 0.58) {
    
    
    
    sig_res <- dplyr::filter(res_tbl, padj < padj.thresh) %>%
        dplyr::arrange(padj)
    
    
    ## Order results by padj values
    top20_sig_genes <- sig_res %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)
    
    ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
    res_table_thres <- res_tbl %>%
        mutate(threshold = padj < padj.thresh & abs(log2FoldChange) >= lfc.thresh) %>% 
        mutate(direction = factor(sign(threshold  * log2FoldChange))) %>%
        mutate(label = ifelse((gene %in% top20_sig_genes) & threshold, gene, NA))
    
    
    ## Volcano plot
    pp.volcano = ggplot(res_table_thres, aes(x = log2FoldChange, 
                                             y = -log10(padj), 
                                             colour = direction,
                                             label = label)) +
        geom_point() +
        geom_text_repel() +
        ggtitle(title) +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        # scale_y_continuous(limits = c(0,50)) +
        scale_color_manual(
            values = c("-1" = muted("blue", 30, 90 ), 
                       "0" = "gray", 
                       "1" = muted("red", 30, 90) ) ) +
        guides (color = "none" ) + 
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25))) + 
        theme_classic(base_size = 13)
    
    return(pp.volcano)
}



read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, function(y) {
        y = as.data.frame(y)
        # y[[1]] = NULL
        return(y)
    })
    names(x) <- sheets
    x
}

write_excel_allsheets <- function(file.list, filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    write.xlsx(file.list, file = filename)
}


add_mouse_geneinfo = function(filename, suffix = "gene_info", entrez.info.file = NULL) {
    
    if (!exists("entrez.gene.condensed")) {
        
        if (is.null(entrez.info.file)) {
            entrez.info.file = here("data/mouse_gene_info/mouse.entrez.info.RDS")
        }
        entrez.gene.all = readRDS(entrez.info.file)
        
        entrez.gene.condensed = do.call(
            rbind,
            lapply(
                entrez.gene.all, 
                function(gene.bunch) {
                    out = sapply(gene.bunch, 
                                 function(x) return(x[c("uid", "name", "description", 
                                                        "summary", "otheraliases", "otherdesignations")])) %>% 
                        t() %>% as.data.frame()
                    return(out)
                } ) ) 
        
        for(cn in colnames(entrez.gene.condensed)) {
            entrez.gene.condensed[[cn]] = unlist(entrez.gene.condensed[[cn]])
        }
    }
    
    outfile = paste0(tools::file_path_sans_ext(filename), "_", suffix, ".xlsx")
    
    file.sheets.list = read_excel_allsheets(filename)
    
    sheet.names = readxl::excel_sheets(filename)
    
    for(sheet.name in sheet.names) {
        file.sheets.list[[sheet.name]] = left_join(
            dplyr::rename(file.sheets.list[[sheet.name]], name = names), 
            entrez.gene.condensed, by = "name")
    }
    
    write_excel_allsheets(file.sheets.list, filename = outfile)
}



read_gene_sets = function(gene.set.path) {
    
    # Open the file connection
    file_conn <- file(gene.set.path, "r")
    
    # Create a list to store the data
    data_list <- list()
    
    # Read the file line by line
    while (TRUE) {
        # Read a line from the file
        line <- readLines(file_conn, n = 1)
        
        # Break the loop if end of file is reached
        if (length(line) == 0) {
            break
        }
        
        # Split the line by tab
        line_elements <- strsplit(line, "\t")[[1]]
        
        # Extract the set name from the first element
        set_name <- line_elements[1]
        
        # Extract the elements belonging to the set
        set_elements <- line_elements[c(-1, -2)]
        
        # Store the set elements in the list
        data_list[[set_name]] <- set_elements
    }
    
    # Close the file connection
    close(file_conn)
    
    # Return the list
    return(data_list)
}

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_mouse_to_human <- function(gene_list){
    
    output = c()
    
    for(gene in gene_list){
        class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
        if(!identical(class_key, integer(0)) ){
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
            for(human_gene in human_genes){
                output = append(output,human_gene)
            }
        }
    }
    
    return (output)
}

