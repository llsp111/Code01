library(shiny)
library(shinythemes)
library(DESeq2)
library(limma)
library(survival)
library(survminer)

ui <- fluidPage(
  theme = shinytheme("sandstone"),
  titlePanel("TCGA表达矩阵分析"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "上传表达矩阵", accept = c(".csv", ".txt")),
      selectInput("gene", "选择基因", choices = NULL),
      actionButton("analyze", "分析"),
      hr(),
      h4("差异基因分析"),
      plotOutput("diffPlot"),
      hr(),
      h4("生存分析"),
      plotOutput("survPlot")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("差异基因分析", dataTableOutput("diffTable")),
        tabPanel("生存分析", dataTableOutput("survTable")),
        tabPanel("示例", 
                 h4("差异基因分析示例"),
                 plotOutput("diffPlotExample"),
                 hr(),
                 h4("生存分析示例"),
                 plotOutput("survPlotExample"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, sep = "\t", row.names = 1)
  })
  
  observe({
    req(data())
    updateSelectInput(session, "gene", choices = colnames(data()))
  })
  
  # 差异基因分析示例数据
  output$diffPlotExample <- renderPlot({
    example_data <- matrix(rnorm(1000), ncol=10)
    colnames(example_data) <- paste0("Sample", 1:10)
    rownames(example_data) <- paste0("Gene", 1:100)
    sample_info <- factor(rep(c("Normal", "Tumor"), each = 5))
    design <- model.matrix(~ sample_info)
    fit <- lmFit(example_data, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=2, number=Inf)
    
    plot(res$logFC, -log10(res$P.Value), 
         xlab="Log2 Fold Change", ylab="-Log10 P-value",
         main="Volcano Plot of Differential Expression (Example)",
         pch=20, col=ifelse(res$adj.P.Val < 0.05, "red", "black"))
  })
  
  # 差异基因分析
  diffAnalysis <- eventReactive(input$analyze, {
    req(data())
    expression_data <- data()
    sample_info <- factor(rep(c("Normal", "Tumor"), each = ncol(expression_data)/2))
    design <- model.matrix(~ sample_info)
    fit <- lmFit(expression_data, design)
    fit <- eBayes(fit)
    topTable(fit, coef=2, number=Inf)
  })
  
  output$diffTable <- renderDataTable({
    req(diffAnalysis())
    as.data.frame(diffAnalysis())
  })
  
  output$diffPlot <- renderPlot({
    req(diffAnalysis())
    res <- diffAnalysis()
    plot(res$logFC, -log10(res$P.Value), 
         xlab="Log2 Fold Change", ylab="-Log10 P-value",
         main="Volcano Plot of Differential Expression",
         pch=20, col=ifelse(res$adj.P.Val < 0.05, "red", "black"))
  })
  
  # 生存分析示例数据
  output$survPlotExample <- renderPlot({
    set.seed(123)
    time <- rexp(100, rate = 0.1)
    status <- sample(0:1, 100, replace = TRUE)
    gene_expr <- rnorm(100)
    median_expr <- median(gene_expr)
    group <- ifelse(gene_expr > median_expr, "High", "Low")
    surv_data <- data.frame(time = time, status = status, group = group)
    surv_fit <- survfit(Surv(time, status) ~ group, data = surv_data)
    
    ggsurvplot(surv_fit, data = surv_data, risk.table = TRUE)
  })
  
  # 生存分析
  survAnalysis <- eventReactive(input$analyze, {
    req(data())
    gene_expr <- data()[, input$gene]
    median_expr <- median(gene_expr)
    group <- ifelse(gene_expr > median_expr, "High", "Low")
    # 使用随机生成的生存时间和状态数据作为示例
    surv_data <- data.frame(time = rnorm(nrow(data()), mean = 1000, sd = 500), 
                            status = sample(0:1, nrow(data()), replace = TRUE),
                            group = group)
    surv_fit <- survfit(Surv(time, status) ~ group, data = surv_data)
    list(surv_fit = surv_fit, surv_data = surv_data)
  })
  
  output$survTable <- renderDataTable({
    req(survAnalysis())
    survAnalysis()$surv_data
  })
  
  output$survPlot <- renderPlot({
    req(survAnalysis())
    ggsurvplot(survAnalysis()$surv_fit, data = survAnalysis()$surv_data, risk.table = TRUE)
  })
}

shinyApp(ui = ui, server = server)


# shiny::runApp('./shiny-app')
