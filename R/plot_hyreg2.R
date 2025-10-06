

#' function to decode which group or observation was classified to which class by the model
#'
#' @description This function can be used to decode the classified classes by the model generates using `hyreg2` or `hyreg2_het` and see,
#'  which group or observation was signed to which class
#'
#' @param data a `dataframe`, which was used to estimate the `model`
#' @param model a flexmix `model`object estimated using [hyreg2()] or [hyreg2_het()]
#' @param id_col `character` string, name of grouping variable, which must be a column of the provided `data`.
#'          the parameter must be specified, if the provided `model` was estimated under control for `groups`
#'
#'
#' @return `dataframe` of two columns, first column named as provided `id_col`  or `"observation"` if `id_col` was not given as
#'          an input. second column named `"mod_comp"` indicating the assigned class for this group or observation
#'
#'
#' @examples
#' # estimate a model using simulated_data_rnorm
#'
#' ### using grouping variable id ####
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 1
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')

#' # use of function give_id
#' give_id(data = simulated_data_norm,
#' model = hyflex_mod,
#' id_col = "id")
#'
#'
#'
#' @author Svenja Elkenkamp & John Grosser
#' @export
#'
#'



give_id <- function(data,
                    model,
                    id_col = NULL) # id must be provided, if model was estimated using a grouping variable
  {

    if(!is.list(model)){
      if(!is.null(id_col)){
        ids_comp <- data.frame(data[,id_col],model@cluster)
        colnames(ids_comp) <- c(id_col,"mod_comp")
      }else{
        data$mod_comp <- model@cluster
        ids_comp <- data.frame(rownames(data),data[,"mod_comp"])
        colnames(ids_comp) <- c("observation","mod_comp")
      }

    }else{
      ids_comp <- model[[length(model)]]
      colnames(ids_comp) <- c(id_col,"mod_comp")
    }

    return(ids_comp)
}





#' plot function to visualize the classification based on the model estimated using `hyreg2` or `hyreg2_het`
#'
#' @description This function can be used to visualize the classification based on the model for different variables.
#'              [ggplot2::ggplot()] is used.
#'
#' @param data a `dataframe`, which was used to estimate the `model` using [hyreg2()] or [hyreg2_het()]
#' @param x `charachter` string, column of `data` to be plotted in x-axis
#' @param y `charachter` string, column of `data` to be plotted in y-axis
#' @param id_col `charachter` sting, grouping variable, same as was given in `model`.
#'            if model was estimated without grouping, see Details
#' @param id_df_model `dataframe` of two columns indicating which group belongs to which class,
#'                  first column named as input `id_col`, second column named `"mod_comp"`.
#'                  this input can be generated using the [give_id()] function, see Details.
#' @param type_to_plot `list` of two `charachter` elements. First: `columnname` of column containing indicator for `type` of `data`,
#'                       Second: `value` of column `type`, that should be used for the plot, see details of [hyreg2()] inputs `type` and `type_cont`,`type_dich`
#' @param colors `charachter` vector, colors to be used in `ggplot`, default `NULL` - than colors are choosen automatically
#'
#'
#'
#' @return `ggplot` object visualizing x against y by classes from the model
#'
#' @details
#' `id_col_df` has to be provided anyway, even if the model was estimated without grouping variable.
#' Since there might be no grouping varibale in the `data`, we recommend to create a new column called `"observation"`
#' in data using the `rownames`/`observationnumbers` as `charachter` values and use this column as
#' input for `id_col` in `plot_hyreg2`, additionally you can then use `id_df_model` =  `give_id(data,model,"observation")`,
#' see example
#'
#' @examples
#' # estimate a model using simulated_data_rnorm
#'
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 2
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "cont",
#'                     id_col = "id"
#')

#'# plotting the variables id against y
#'plot_hyreg2(data = simulated_data_norm,
#'           x = "id",
#'           y = "y",
#'           id_col = "id",
#'           id_df_model = give_id(data = simulated_data_norm,
#'                                 model = hyflex_mod,
#'                                 id = "id"))
#'

#'
#' @author Svenja Elkenkamp & John Grosser
#' @importFrom ggplot2 ggplot
#' @export
#'
#'

# include colot option
plot_hyreg2 <- function(data,
                       x,
                       y,
                       id_col,
                       id_df_model,  # you can use give_id() to generate id_df
                       type_to_plot = NULL, #list of two elements list("type","TTO")
                       colors = NULL # optional colour vector
){


  colnames(id_df_model) <- c(id_col,"mod_comp")
  data <- merge(data, id_df_model, by = id_col)

  if(!is.null(type_to_plot)){
    data <- data[(data[,type_to_plot[[1]]]) == type_to_plot[[2]],]
  }

# https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html
#require("ggplot2")

 p <-  ggplot2::ggplot(mapping = ggplot2::aes(x = data[,x], y = data[,y], color = as.character(data$mod_comp)))
 p <- p + ggplot2::geom_point(size = 1, alpha = 0.7)
 p <- p +  ggplot2::labs( title = "classification",
          x = paste0(x),
          y = paste0(y),
          color = "Class")
 p <- p + ggplot2::geom_jitter(width = 0.1, height = 0, size = 1) #alpha = 0.2) # verschieben der Punkte leicht zur Seite
 p <- p + ggplot2::theme_minimal()

  # colours
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

 return(p)

}







