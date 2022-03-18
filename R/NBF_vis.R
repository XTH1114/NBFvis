#' Neighborhood-Based Featurization Visualization
#'
#' The function could apply UMAP to the featurization matrix of a special omics dataset. \cr
#' This matrix could be a neighborhood matrix, which is the combination of the quantile matrix and network matrix (from quantile_matrix and network_matrix, respectively) \cr
#' Then, a shiny app will visualize the UMAP embeddings of the spatial omics data.
#'
#' @param matrix a statistics matrix. It should be rescaled. \cr
#' For example, it could be the rescaled quantiles network centralities matrix.
#' @param origin_data a data frame. It is the origin dataset of cell expressions, like proteins and genes. \cr
#' Column \strong{index}, \strong{Group}, \strong{x_center} and \strong{y_center} should be included.
#' @param var_names the names of columns about expressions like proteins or genes in the argument \strong{origin_data}.
#' @return an interactive shiny app for the UMAP embeddings.
#' @importFrom umap umap
#' @importFrom dplyr mutate select filter everything left_join group_by
#' @importFrom tidyr spread
#' @importFrom magrittr %>% %<>%
#' @importFrom heatmaply heatmaply
#' @rawNamespace import(ggplot2, except = last_plot)
#' @import shiny shinythemes plotly RColorBrewer
#' @export
#' @examples
#' # it may take a few seconds to generate the Shiny App.
#' #library(dplyr)
#' #data(patient4)
#' #data(centrality_patient4)
#' #data(Quantiles_patient4)
#' #neighborhood_info_patient4 <- cbind(
#' #Quantiles_patient4 %>% select(-index, -n_neighbor),
#' #scale(centrality_patient4 %>% select(-index)))
#'
#' #NBF_vis(
#'  # matrix = neighborhood_info_patient4,
#'   #origin_data = patient4,
#'   #var_names = colnames(patient4)[17:57])

NBF_vis <- function (matrix, origin_data, var_names){
  complete_index <- complete.cases(matrix)
  matrix <- matrix[complete_index, ]
  origin_data <- origin_data[complete_index, ]
  all_umap <- umap(matrix)
  all_layout <- all_umap$layout
  colnames(all_layout) <- c("x", "y")
  all_umap_plot <- all_layout %>% data.frame() %>% mutate(Group = factor(origin_data$Group),
                                                          index = origin_data$index)
  centers <- origin_data[, c("index", "x_center", "y_center")]
  km_embedding <- left_join(centers, all_umap_plot)
  embedding <- km_embedding %>% group_by(Group)
  embedding1 <- highlight_key(embedding, ~index)
  scat3d <- plot_ly(embedding1, hoverinfo = "text", legendgroup = ~Group) %>%
    add_markers(x = ~x, y = ~y, color = ~Group, text = ~Group, marker = list(size = 4)) %>%
    layout(yaxis = list(title = ""), xaxis = list(title = ""))
  scat <- plot_ly(embedding1, hoverinfo = "text", legendgroup = ~Group,
                  showlegend = F) %>% add_markers(x = ~x_center, y = ~y_center,
                                                  color = ~Group, text = ~Group,
                                                  marker = list(size = 4)) %>%
    layout(title = "UMAP Embedding Plot and Spatial Plot",
           yaxis = list(title = ""), xaxis = list(title = ""))
  origin_data1_plotly <- subplot(scat3d, scat, titleX = TRUE) %>%
    highlight(on = "plotly_selected", off = "plotly_deselect")
  ui <- fluidPage(theme = shinytheme("cerulean"), titlePanel("Neighborhood-Based Featurization",
                                                             "NBFvis Shiny"), verbatimTextOutput("debug"),mainPanel(plotlyOutput("scatters",
                                                                                                                                 height = "600px", width = "1200px"), sliderInput("clusters",
                                                                                                                                                                                  "The number of clusters", value = 5, min = 1, max = 12),
                                                                                                                    actionButton("action", "Confirm"), plotlyOutput("scatters_2",
                                                                                                                                                                    height = "600px", width = "1200px"), uiOutput("selecting_clusters"),
                                                                                                                    uiOutput("selecting_groups"), tabsetPanel(tabPanel("Heatmap",
                                                                                                                                                                       numericInput("num_variables", "The number of the most significant variables shown in the heatmap",
                                                                                                                                                                                    value = 10, min = 1, max = length(var_names)),
                                                                                                                                                                       plotlyOutput("heatmap1", width = "1200px")), tabPanel("Structure map",
                                                                                                                                                                                                                             plotlyOutput("structure", width = "1200px")), tabPanel("histograms",
                                                                                                                                                                                                                                                                                    selectInput("selected_variable", "Select a varible: ",
                                                                                                                                                                                                                                                                                                choices = var_names, selected = var_names[1]),
                                                                                                                                                                                                                                                                                    plotlyOutput("hist", width = "1200px"))), plotlyOutput("scatters_3",
                                                                                                                                                                                                                                                                                                                                           height = "600px", width = "1200px")))
  server <- function(input, output, session) {
    cluster_index <- eventReactive(input$action, {
      kmeans(km_embedding %>% select(x, y), centers = input$clusters,
             nstart = 30) %>% .$cluster
    })
    filtered_data <- reactive({
      norm_data <- origin_data[origin_data$Group %in% as.character(input$selected_groups),
      ] %>% select(var_names) %>% scale()
      location_data <- origin_data[origin_data$Group %in%
                                     as.character(input$selected_groups), ] %>% select(index)
      temp_data <- cbind(index = location_data, norm_data)
      temp_data$clusters <- cluster_index()[origin_data$Group %in%
                                              as.character(input$selected_groups)]
      temp_data
    })
    output$scatters <- renderPlotly({
      origin_data1_plotly
    })
    output$selecting_groups <- renderUI({
      checkboxGroupInput("selected_groups", label = "Select Groups",
                         choices = levels(factor(origin_data$Group)),
                         inline = T)
    })
    output$scatters_2 <- renderPlotly({
      validate(
        need(cluster_index(), "Please press the Confirm button.")
      )
      embedding2 <- highlight_key(km_embedding %>% mutate(cluster = factor(cluster_index())) %>%
                                    group_by(cluster), ~index)
      scat3d_2 <- plot_ly(embedding2, hoverinfo = c("text"),
                          legendgroup = ~cluster) %>%
        add_markers(x = ~x,y = ~y, color = ~cluster,
                    text = ~cluster,marker = list(size = 4),
                    colors = brewer.pal(input$clusters,"Set2")) %>%
        layout(title = "K-means Clustering Plot and Spatial Plot",
               yaxis = list(title = ""), xaxis = list(title = ""))
      scat_2 <- plot_ly(embedding2, hoverinfo = c("text"),
                        legendgroup = ~cluster, showlegend = F) %>% add_markers(x = ~x_center,
                                                                                y = ~y_center,
                                                                                color = ~cluster, marker = list( size = 4),
                                                                                colors = brewer.pal(input$clusters,
                                                                                                    "Set2")) %>%
        layout(yaxis = list(title = ""), xaxis = list(title = ""))
      origin_data2_plotly <- subplot(scat3d_2, scat_2,
                                     titleX = TRUE) %>% highlight(on = "plotly_selected",
                                                                  off = "plotly_deselect")
      origin_data2_plotly
    })
    output$selecting_clusters <- renderUI({
      checkboxGroupInput("selected_clusters", label = "Select K-means Clusters",
                         choices = 1:input$clusters, inline = T)
    })

    #select K-means clusters
    output$selecting_clusters <- renderUI({
      checkboxGroupInput("selected_clusters", label = "Select K-means Clusters",
                         choices = 1:input$clusters, inline = T)
    })
    #heatmap
    output$heatmap1 <- renderPlotly({
      validate(
        need(cluster_index(), "Please press the Confirm button.")
      )
      temp_data <- filtered_data()
      important_variables <- aggregate(. ~ clusters, select(temp_data,
                                                            -index), FUN = median, na.action = na.pass) %>%
        filter(clusters %in% as.character(input$selected_clusters)) %>%
        select(-clusters) %>% sapply(FUN = var, na.rm = T) %>%
        sort(decreasing = T) %>% names() %>% .[1:input$num_variables]
      temp_data2 <- temp_data %>% filter(clusters %in%
                                           as.character(input$selected_clusters))
      heatmap_data <- aggregate(temp_data2 %>% select(one_of(important_variables)),
                                by = list(group = temp_data2$clusters), median,
                                na.rm = T)
      rownames(heatmap_data) <- as.character(heatmap_data$group)
      heatmap_data %<>% select(-group)
      heatmaply(heatmap_data[, order(names(heatmap_data))],
                column_text_angle = 90, Rowv = F, Colv = T, hclust_method = NULL, labRow = paste0("K-means Cluster",input$selected_clusters))
    })
    #histogram
    output$hist <- renderPlotly({
      validate(
        need(cluster_index(), "Please press the Confirm button.")
      )
      temp_data <- filtered_data()
      temp_data$clusters %<>% factor
      temp_data <- temp_data %>% filter(clusters %in%
                                          input$selected_clusters)
      ggplotly(ggplot(temp_data) + geom_histogram(aes_string(input$selected_variable,
                                                             y = "..density..", fill = "clusters"), alpha = 0.6,
                                                  position = "identity") + scale_fill_manual(values = brewer.pal(input$clusters,
                                                                                                                 "Set2")[as.numeric(input$selected_clusters)])+
                 guides(fill=guide_legend(title="Cluster")))
    })
    #structure plot
    output$structure <- renderPlotly({
      validate(
        need(cluster_index(), "Please press the Confirm button.")
      )
      temp_data <- filtered_data()
      temp_data <- temp_data %>%
        mutate(Group = origin_data[origin_data$Group %in% as.character(input$selected_groups),"Group"])
      temp_data$clusters %<>% factor
      temp_data <- temp_data %>% filter(clusters %in%
                                          input$selected_clusters)
      temp_data$Group <- factor(temp_data$Group)
      temp_data$clusters <- factor(temp_data$clusters %>% as.character)

      hclust_group_cluster <- temp_data %>%
        group_by(clusters) %>%
        mutate(clusters_num = n()) %>%
        group_by(clusters,Group) %>%
        summarise(Group_num = n(),clusters_num = mean(clusters_num)) %>%
        mutate(ratio = Group_num/clusters_num) %>%
        as.data.frame() %>%
        select(clusters,Group,ratio) %>% spread("Group","ratio")
      hclust_group_cluster_order <-  hclust_group_cluster %>%
        .[,-1] %>%
        dist() %>%
        hclust %>%
        .$order


      temp_data$clusters <- factor(temp_data$clusters,
                                   levels = as.character(hclust_group_cluster$clusters)[hclust_group_cluster_order],
                                   ordered = T)



      ggplotly(ggplot(temp_data)+
                 geom_bar(aes(clusters,fill = Group),size = 0.5,width = 1, position = "fill")+
                 coord_flip()+
                 guides(fill=guide_legend(title="")) +
                 scale_fill_brewer(palette = "Set2"))

    })

    output$debug <- renderPrint({
      temp_data <- filtered_data()
      temp_data <- temp_data %>%
        mutate(Group = origin_data[origin_data$Group %in% as.character(input$selected_groups),"Group"])
      temp_data$clusters %<>% factor
      temp_data <- temp_data %>% filter(clusters %in%
                                          input$selected_clusters)
      temp_data$Group <- factor(temp_data$Group)
      temp_data$clusters <- factor(temp_data$clusters %>% as.character)

      hclust_group_cluster <- temp_data %>%
        group_by(clusters) %>%
        mutate(clusters_num = n()) %>%
        group_by(clusters,Group) %>%
        summarise(Group_num = n(),clusters_num = mean(clusters_num)) %>%
        mutate(ratio = Group_num/clusters_num) %>%
        as.data.frame() %>%
        select(clusters,Group,ratio) %>% spread("Group","ratio")
      hclust_group_cluster_order <-  hclust_group_cluster %>%
        .[,-1] %>%
        dist() %>%
        hclust %>%
        .$order
      hclust_group_cluster
    })


    output$scatters_3 <- renderPlotly({
      embedding3 <- highlight_key(km_embedding %>% mutate(cluster = factor(cluster_index()),
                                                          !!input$selected_variable := origin_data[, input$selected_variable]) %>%
                                    filter(Group %in% input$selected_groups) %>%
                                    group_by(cluster), ~index)
      scat3d_3 <- plot_ly(embedding3, hoverinfo = c("text"),
                          legendgroup = ~cluster) %>% add_markers(x = ~x,
                                                                  y = ~y, color = as.formula(paste0('~', input$selected_variable)), text = ~cluster,marker = list(size = 4)) %>%
        layout(title = paste0("UMAP embeddings and expression plot of ", input$selected_variable),
               yaxis = list(title = ""), xaxis = list(title = ""), legend=list(title=list(text=input$selected_variable)))
      scat_3 <- plot_ly(embedding3, hoverinfo = c("text"),
                        legendgroup = ~cluster) %>% add_markers(x = ~x_center, y = ~y_center, color = as.formula(paste0('~', input$selected_variable)), marker=list(size = 4)) %>%
        layout(showlegend = F,yaxis = list(title = ""), xaxis = list(title = ""))
      origin_data3_plotly <- subplot(scat3d_3, scat_3,
                                     titleX = TRUE) %>% highlight(on = "plotly_selected",
                                                                  off = "plotly_deselect")
      origin_data3_plotly
    })

  }
  shinyApp(ui, server)
}

