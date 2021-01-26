library(shiny)
library(tidyverse)
library(plotly)
library("vroom")
library(shinybusy)

####################

# library(rsconnect)
# options(repos = BiocManager::repositories())
# rsconnect::deployApp(".")


###### All stages
vroom(file = "data2/Expression.tsv",delim = "\t",) -> expression
factor(expression$Stage, levels= c("Miracidia","Sporocysts","Cercariae","Somula","Males","Females")) -> expression$Stage

read.table("data2/CV.tsv",sep = "\t", header = T) -> CV_wide
gather(CV_wide,"Normalization","CV",-gene_ID) -> CV_long

CV_wide[order(CV_wide$DESeq2, decreasing = F),] -> CV_wide

CV_wide$gene_ID -> genes
CV_wide[,c(1,3,2,4,5)] -> CV_wide

###### Filtered stages
vroom(file = "data1/Expression2.tsv",delim = "\t",) -> expression2
factor(expression2$Stage, levels= c("Cercariae","Somula","Males","Females")) -> expression2$Stage

read.table("data1/CV2.tsv",sep = "\t", header = T) -> CV_wide2
gather(CV_wide2,"Normalization","CV",-gene_ID) -> CV_long2

CV_wide2[order(CV_wide2$DESeq2, decreasing = F),] -> CV_wide2

CV_wide2$gene_ID -> genes2
CV_wide2[,c(1,3,2,4,5)] -> CV_wide2


ui <- navbarPage("App S. mansoni",
    # Application title
    tabPanel("All stages",
    
    fluidRow(
        column(12,align="center",
               selectInput(inputId = "select_gene", 
                                                   label = "Select gene ID to plot",
                                                   choices = genes,
                                                   multiple = F,
                                                   selected = genes[1]))),
    fluidRow(
        column(3,align="left",
            DT::dataTableOutput("cv_df"),
            downloadButton("downloadData", "Download")
        ),
        column(8,offset = 1,align="right",
             plotlyOutput("lineplot",height = 300),
             # downloadLink("downloadPlot", "Download Plot"),
             plotlyOutput("cv_plot",height = 200),
             ))
    ),
    
    tabPanel("Filtered stages",
             
             fluidRow(
                 column(12,align="center",
                        selectInput(inputId = "select_gene2", 
                                    label = "Select gene ID to plot",
                                    choices = genes2,
                                    multiple = F,
                                    selected = genes2[1]))),
             fluidRow(
                 column(3,align="left",
                        DT::dataTableOutput("cv_df2"),
                        downloadButton("downloadData2", "Download")
                 ),
                 column(8,offset = 1,align="right",
                        plotlyOutput("lineplot2",height = 300),
                        # downloadLink("downloadPlot2", "Download Plot"),
                        plotlyOutput("cv_plot2",height = 200),
                 ))
    )
    )



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    ###### All 
    output$lineplot <- renderPlotly({
        show_modal_spinner()
        req(input$select_gene)
        expression[expression$gene_ID == input$select_gene,] -> selected_df
        remove_modal_spinner()
        
        ggplot(selected_df, aes(x = Stage, y = expression, color=method , sample= Sample)) +
            facet_wrap(~method,scales = "free",nrow = 1)+
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            # stat_summary(fun=median, geom="point", size=1, color="red", show.legend = F)+
            stat_summary(fun=median, colour="red", geom="line", aes(group = 1))+
            theme(plot.title = element_text(face = "bold",hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.text.x = element_text(angle = 45,hjust=1),
                  legend.position = "none")+
            ggtitle(input$select_gene)+
            xlab("")+
            ylab("Expression")
    })
    
    output$cv_plot <- renderPlotly({
        req(input$select_gene)
        CV_long[CV_long$gene_ID %in% input$select_gene,] -> CV_long_sel
        
        ggplot(CV_long_sel, aes(x = Normalization,y = CV,color=Normalization, gene= gene_ID)) +
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            theme(plot.title = element_text(face = "bold",hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.text.x = element_text(angle = 45,hjust=1),
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

    ############ Filtered
    
    output$lineplot2 <- renderPlotly({
        show_modal_spinner()
        req(input$select_gene2)
        expression2[expression2$gene_ID == input$select_gene2,] -> selected_df
        remove_modal_spinner()
        
        ggplot(selected_df, aes(x = Stage, y = expression, color=method , sample= Sample)) +
            facet_wrap(~method,scales = "free",nrow = 1)+
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            # stat_summary(fun=median, geom="point", size=1, color="red", show.legend = F)+
            stat_summary(fun=median, colour="red", geom="line", aes(group = 1))+
            theme(plot.title = element_text(face = "bold",hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.text.x = element_text(angle = 45,hjust=1),
                  legend.position = "none")+
            ggtitle(input$select_gene2)+
            xlab("")+
            ylab("Expression")
    })
    
    output$cv_plot2 <- renderPlotly({
        req(input$select_gene2)
        CV_long2[CV_long2$gene_ID %in% input$select_gene2,] -> CV_long_sel
        
        ggplot(CV_long_sel, aes(x = Normalization,y = CV,color=Normalization, gene= gene_ID)) +
            geom_point(alpha=1, size=1) +
            scale_color_brewer(palette="Dark2")+
            theme(plot.title = element_text(face = "bold",hjust = 0.5),
                  strip.text.x = element_text(size = 12,face = "bold"),
                  axis.text.x = element_text(angle = 45,hjust=1),
                  legend.position = "none")+
            ggtitle("")+
            xlab("")+
            ylab("CV")
        
    })
    
    output$cv_df2 <- DT::renderDataTable(
        (CV_wide2 %>% mutate(across(is.numeric, ~ round(., 4)))) ,rownames= FALSE
    ) 
    
    
    #Download the PUL file
    output$downloadData2 <- downloadHandler(
        filename = function() {
            paste("CV2.tsv", sep = "")
        },
        content = function(file) {
            write.table(CV_wide2, file, row.names = FALSE, sep = "\t")
        })
    
}





# Run the application 
shinyApp(ui = ui, server = server)
