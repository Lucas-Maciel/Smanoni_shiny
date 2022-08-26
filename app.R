library(shiny)
library(tidyverse)
library(plotly)
library("vroom")
library(shinybusy)
library(dplyr)

####################

###### All stages
vroom(file = "data2/Expression.tsv",delim = "\t",) -> expression
expression$stages2 <- expression$Stage
unique(expression$method)
factor(expression$method, levels= c("DESeq2","TMM.CPM","UQ","TPM")) -> expression$method

#Grouping and chaning names
str_replace_all(expression$stages2,"Miracidia","M/S") -> expression$stages2
str_replace_all(expression$stages2,"Sporocysts","M/S") -> expression$stages2
str_replace_all(expression$stages2,"Cercariae","C") -> expression$stages2
str_replace_all(expression$stages2,"Somula","So") -> expression$stages2
str_replace_all(expression$stages2,"Males","M") -> expression$stages2
str_replace_all(expression$stages2,"Females","F") -> expression$stages2

factor(expression$Stage, levels= c("Miracidia","Sporocysts","Cercariae","Somula","Males","Females")) -> expression$Stage
factor(expression$stages2, levels= c("M/S","C","So","M","F")) -> expression$stages2

read.table("data2/CV.tsv",sep = "\t", header = T) -> CV_wide
gather(CV_wide,"Normalization","CV",-gene_ID) -> CV_long

CV_wide[order(CV_wide$DESeq2, decreasing = F),] -> CV_wide

CV_wide$gene_ID -> genes
CV_wide[,c(1,3,2,5,4)] -> CV_wide

###### Filtered stages
# vroom(file = "data1/Expression2.tsv",delim = "\t",) -> expression2
# factor(expression2$Stage, levels= c("Cercariae","Somula","Males","Females")) -> expression2$Stage
# 
# read.table("data1/CV2.tsv",sep = "\t", header = T) -> CV_wide2
# gather(CV_wide2,"Normalization","CV",-gene_ID) -> CV_long2
# 
# CV_wide2[order(CV_wide2$DESeq2, decreasing = F),] -> CV_wide2
# 
# CV_wide2$gene_ID -> genes2
# CV_wide2[,c(1,3,2,4,5)] -> CV_wide2

url <- a("https://doi.org/10.1038/s41598-021-96055-7", href="https://doi.org/10.1038/s41598-021-96055-7")

ui <- navbarPage("Schistosoma mansoni reference genes",
    # Application title
    tabPanel("VerjoLab",
             
    fluidRow(column(12,align="left",
                    h5("This site shows the expression of Schistosoma mansoni protein-coding genes (Smps) and long non-coding RNAs (SmLINC, SmLNCA, SmLNCS RNAs) that comprise the dataset which was used to determine stable reference genes (Smps) across six life-cycle stages of the parasite (Silveira et al., Assessment of reference genes at six different developmental stages of Schistosoma mansoni for quantitative RT-PCR, Scientific Reports 2021 ",
                       url,").")
        
    )),
    fluidRow(column(12,align="left",
                    h5("As explained in the Silveira et al. 2021 paper, after reads counting with RSEM we performed a minimal pre-filtering to keep only genes that had at least 10 reads total when adding all stages. This resulted in 13,624 protein-coding genes (out of 14,548 genes in the v7.1 transcriptome (WBPS14), PMID: 27899279) and in 9,391 lncRNAs (out of 16,583 lncRNAs, PMID: 31572441) that were considered for further analyses. These genes are shown in the present site. Four different methods were used to normalize the expression across all samples from the six different life-cycle stages, namely DESeq2, Trimmed mean of M-values (TMM.CPM), Upper Quartile (UQ), and Transcripts per Million (TPM).")
    )),
    fluidRow(column(12,align="left",
                    h5("You can search for your gene of interest (Smp or lncRNA) by typing its gene ID name below.")
    )),
             
    fluidRow(
        column(12,align="center",
               selectInput(inputId = "select_gene", 
                                                   label = "Type gene ID to plot",
                                                   choices = genes,
                                                   multiple = F,
                                                   selected = genes[1]))),
    fluidRow(
        column(7,align="center",height = 1000,#offset = 1
               plotlyOutput("lineplot",height = 300) ,
               plotlyOutput("lineplot_teste",height = 320) 
               # downloadLink("downloadPlot", "Download Plot"),
        ),
        column(4,align="center",
               plotlyOutput("cv_plot",height = 170),
               tags$h4("Coeficient of Variation (CV) for each method"),
               DT::dataTableOutput("cv_df"),
               downloadButton("downloadData", "Download")
        ))
    )
    )

#

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    ###### All 
    
    dataInput <- reactive({
        req(input$select_gene)
        show_modal_spinner()
        expression[expression$gene_ID == input$select_gene,] -> selected_df
        remove_modal_spinner()
        selected_df
        })
    
    output$lineplot <- renderPlotly({
        # show_modal_spinner()
        req(input$select_gene)
        # remove_modal_spinner()
        dataInput() -> selected_df
        selected_df[selected_df$method == "DESeq2" | selected_df$method == "TMM.CPM", ] -> selected_df 
        g1 <- ggplot(selected_df, 
                     aes(x = stages2, y = expression, color=Stage , sample= Sample)) +
            facet_wrap(~method,scales = "free",nrow = 1)+
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            # stat_summary(fun=median, geom="point", size=1, color="red", show.legend = F)+
            stat_summary(fun=median, colour="red", geom="line", aes(group = 1))+
            theme(plot.title = element_text(face = "bold",vjust = 3,hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(face="bold",size = 11),
                  axis.text.y = element_text(face="bold",size=11),
                  legend.position = "none")+
            ggtitle(input$select_gene)
            
        y <- list(
            title = "Expression",title_y=-1
        )
        ggplotly(g1) %>%  layout(title = list(y = .99), 
                                 yaxis =  y,
                                 margin = list(l = 75)) #%>%
            #layout(legend = list(orientation = 'h',y=-0.1))
    })
    
    output$lineplot_teste <- renderPlotly({
        req(input$select_gene)
        dataInput() -> selected_df
        selected_df[selected_df$method == "UQ" | selected_df$method == "TPM" ,] -> selected_df
        
        g2 <- ggplot(selected_df, aes(x = stages2, y = expression, color=Stage , sample= Sample)) +
            facet_wrap(~method,scales = "free",nrow = 1)+
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            stat_summary(fun=median, colour="red", geom="line", aes(group = 1))+
            theme(plot.title = element_text(face = "bold",vjust = 3,hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(face="bold",size = 11),
                  axis.text.y = element_text(face="bold",size=11),
                  legend.position = "bottom")
        
        y <- list(
            title = "Expression",title_y=-1
        )
        ggplotly(g2) %>%  layout( yaxis =  y,
                                 margin = list(l = 75)) %>%
            layout(legend = list(orientation = 'h',y=-0.1))
    })
    
    
    output$cv_plot <- renderPlotly({
        req(input$select_gene)
        CV_long[CV_long$gene_ID %in% input$select_gene,] -> CV_long_sel
        
        ggplot(CV_long_sel, aes(x = Normalization,y = CV,color=Normalization, gene= gene_ID)) +
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            theme(plot.title = element_text(face = "bold",hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.text.x = element_text(angle = 0,hjust=1,face="bold"),
                  legend.position = "none")+
            ggtitle("")+
            xlab("")+
            ylab("CV")
        
    })
    
    output$cv_df <- DT::renderDataTable(
        (CV_wide %>% mutate(across(is.numeric, ~ round(., 4)))) ,rownames= FALSE
    ) 
    
    #Download the PUL file
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("CV.tsv", sep = "")
        },
        content = function(file) {
            write.table(CV_wide, file, row.names = FALSE, sep = "\t")
        })
}





# Run the application 
shinyApp(ui = ui, server = server)
