#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(shinyBS)



# Define UI for application that draws a histogram
ui <- fluidPage(setBackgroundColor("ivory"),
  headerPanel('Five-year survival probability for incident COPD patients'),br(),
    h4('This web calculator is intended to be used ONLY by HEALTHCARE PROFESSIONALS and relates ONLY to INCIDENT COPD (i.e. COPD diagnosed very recently), we do not recommend its use by patients.  This model has NOT BEEN VALIDATED USING DATA FROM COUNTRIES OTHER THAN ENGLAND. Further details on the model and its development provided below.'),
      h5('Full details of the development of this model are available at LINK. In short, this model was trained using data from English GP records from the Vision system made available through CPRD-GOLD. It was trained using data on 47,964 patients and tested on data from 12,096 patients.'),br(),
  fluidRow(
      column(4,
      numericInput('age', 'Age in years', 67.7),
    selectInput('gender', 'Gender',c('Female','Male')),
    selectInput('smoke', 'Smoking status', c('Current smoker','Ex-smoker','Never smoker')),
    checkboxInput("bmi_NR", "BMI not recorded (i.e. not known)"),
    conditionalPanel( condition = "input.bmi_NR == false",
                    numericInput('bmi', 'BMI', 26)),
    checkboxInput("fev1_NR", "FEV1 % predicted not recorded  (i.e. not known)"),
    conditionalPanel( condition = "input.fev1_NR == false",
                    numericInput('fev1', 'FEV1 % predicted', 64.6))
      ),
    column(4,
           tags$h4("Patient conditions (please select ALL that apply)"),br(),
    checkboxInput("ap", "Alcohol abuse (ever)"),
    bsTooltip(id='ap',title='Has the patient ever had a diagnosis of health problems related to excessive alcohol use?',placement = 'left',trigger='hover'),
    checkboxInput("ast", "Asthma (for at least two years, medication this year)"),
        bsTooltip(id='ast',title='Has the patient received a diagnosis of asthma prior to the two years before their COPD diagnosis that has required medication in the last 12 months?',placement = 'left',trigger='hover'),
    checkboxInput("atr", "Atrial fibrillation (ever)"),
    checkboxInput("can", "Cancer (last five years)"),
        bsTooltip(id='can',title='Has the patient had a cancer diagnosis recorded in the last 5 years?',placement = 'left',trigger='hover'),
    checkboxInput("rhe", "Connective tissue disorder (ever, e.g. rheumatoid arthritis)"),
 checkboxInput("con", "Constipation"),
        bsTooltip(id='con',title='Has this patient had 4 or more laxative prescriptions in the last 12 months?',placement = 'left',trigger='hover'),
    checkboxInput("dep", "Depression"),
        bsTooltip(id='dep',title='Has this patient had a depression code recorded in the last 12 months? OR 4 or more anti-depressant prescriptions (excluding low dose tricyclics) in last year?',placement = 'left',trigger='hover')
 ),
    column(4,    checkboxInput("dia", "Diabetes (ever)"),
        bsTooltip(id='dia',title='Has the patient ever been diagnosed as diabetic (Type 1 or 2)?',placement = 'left',trigger='hover'),
    checkboxInput("epi", "Epilepsy (treated in last 12 months)"),
    bsTooltip(id='epi',title='	Has the patient ever been diagnosed as having epilepsy AND has received at least 1 antiepileptic prescription in the last 12 months?',placement = 'left',trigger='hover'),
    checkboxInput("hf", "Heart failure (ever)"),
    checkboxInput("ibd", "Inflammatory bowel disease (ever)"),
    checkboxInput("ibs", "Irritable bowel syndrome (ever)"),
    bsTooltip(id='ibs',title='Has this patient ever had irritable bowel syndrome? OR have they received 4 or more antispasmodic prescriptions in the last 12 months?',placement = 'left',trigger='hover'),
    checkboxInput("pvd", "Peripheral vascular disorder (ever)"),
    checkboxInput("psm", "Psychoactive substance misuse (ever - not alcohol)"),
    bsTooltip(id='psm',title='Has this patient ever had a problem with psychoactive substance misuse?',placement = 'left',trigger='hover'),
    checkboxInput("scz", "Scizophrenia / bipolar (ever)"),
    bsTooltip(id='scz',title='Has this patient ever been diagnosed with either schizophrenia or bipolar? Or have they ever been prescribed lithium?',placement = 'left',trigger='hover'),
    checkboxInput("str", "Stroke or Transient ischaemic attack (ever)"),br(),actionButton(inputId = "calculate", label = "Calculate survival probability")
    )
  ),hr(),
  mainPanel(
    strong(span(textOutput("note"), style="color:red"))
  )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
isValid_fev1<-reactive({
  !is.null(input$fev1) && input$fev1>=0 && input$fev1<=100
})
    
isValid_age<-reactive({
  !is.null(input$age) && input$age>34 && input$age<106
})

vals <- reactiveValues()

output$note <- renderText({
    
    if(!isValid_age()){'Age needs to be between 35 - 105 years.'} else {
        
        if(input$calculate>=1){
        
            logit <- 1.874370047 + (-0.085820725)*(input$age-67.7) + (-0.001089517)*((input$age-67.7)^2) +
             (0.036567533)*(input$bmi-26) + (-0.002428287)*((input$bmi-26)^2) + (0.013684613)*(input$fev1-64.6) +
            (-0.000277955)*((input$fev1-64.6)^2) +(0.264175108)*(input$gender=='Female') +
            (0.641164355)*(input$smoke=='Never smoker') + (0.402440002)*(input$smoke=='Ex-smoker') +
            (-0.504466266)*input$bmi_NR + (-0.540682601)*input$fev1_NR + (-0.813714263)*input$can +
            (-0.811468347)*input$hf + (-0.720885687)*input$ap + (-0.522215569)*input$pvd +
            (-0.432218107)*input$atr + (0.232582125)*input$ast + (-0.390370747)*input$dia +
            (-0.351968752) *input$str + (-0.379362972) * input$epi + (0.292636023) *input$ibs +
            (-0.412296068) *input$scz + (-0.317463786) *input$psm + (-0.265865416) *input$ibd +
            (-0.286588747) *input$con + (-0.291772449) *input$dep + (-0.245848919) *input$rhe
            
            paste('Estimated probability of five-year survival is: ',round((exp(logit)/(1+exp(logit)))*100),' %')
        
        } else {''}
    }
    
})


}



# Run the application 
shinyApp(ui = ui, server = server)
