\name{cfw.pheno}
\docType{data}
\alias{cfw.pheno}
\title{Phenotype data for 1,092 CFW male mice}
\description{
  The data frame contains physiological and behavioural phenotype data
  on 1,092 mice from the Carworth Farms White (CFW) outbred mouse stock.
  All mice are males. 
}
\usage{data(cfw.pheno)}
\format{
  Individual columns of the table \code{cfw.pheno} are as follows:
  
  \describe{
    \item{id}{Unique number assigned to each mouse.}

    \item{round}{Refers to a shipment of mice. There are at most 48 mice
      in a single shipment.}

    \item{cageid}{Each mouse is placed in a cage with 3 other mice.}

    \item{FCbox, PPIbox, methcage}{Boxes and cages used for behavioral
      testing.}

    \item{glucoseage, methage, FCage, PPIage, sacage}{Age of mouse (in
      days) when glucose levels are measured, when various behavioral
      tests are started, and when mouse is sacrificed.}
    
    \item{bw0, bw1, bw2, bw3, PPIweight, sacweight}{Body weights
      recorded at various time points. \code{bw0} is the body weight
      measured upon arrival in Chicago, when tail is
      snipped. \code{bw1}, \code{bw2} and \code{bw3} are body weights
      taken during methamphetamine sensivity testing. \code{PPIweight}
      is the body weight measured during prepulse inhibition
      testing. \code{sacweight} is the body weight measured when the
      mouse is sacrificed.}

    \item{BMD}{Bone-mineral density.}

    \item{TA, EDL, gastroc, plantaris, soleus}{Measured weights of
       various muscles(in mg).}
    
    \item{tibia}{Length of the tibia bone (in mm).}
     
    \item{abnormalbone}{Refers to "bone health", in which 0 = healthy
      looking bone, and 1 = abnormal bone, white and swollen, which,
      might be expression of osteopetrosis. Note that this is not the
      same "abnormal bone" phenotype that we used for QTL mapping.}
  
    \item{experimenters}{Initials of the people who handled the mice.}

    \item{fastglucose}{Glucose levels in plasma after fasting. Measured
      using glucometer once at dissection, after fasting overnight.
      Units are milligrams per deciliter (mg/dL), which is the SI unit
      for measuring glucose in blood.}

    \item{PreTrainD1}{Average proportion of freezing on day 1 during the
      pre-training interval (30-180 seconds) before exposure to tones
      and shocks ("pre-training freezing").}
    
    \item{AvToneD1}{Average proportion of freezing on the first day
      during exposure to the conditioned stimulus
      ("freezing to tone on day 1").}
      
    \item{AvContextD2}{Average proportion of time freezing in the 30–180
      second interval on the second day in conditions identical to the
      first day ("freezing to same context").}
    
    \item{AvAltContextD3}{Average proportion of time freezing over the
      30-180 interval on the third day in an altered setting.}
      
    \item{AvToneD3}{Average proportion of time freezing on the third day
      in the altered setting during the 30-second intervals in which the
      tones are presented ("freezing to cue").}

    \item{pp3PPIavg, pp6PPIavg, pp12PPIavg}{Average of the
      "inhibition intensity" taken as the ratio of the prepulse response
      during the 3, 6 and 12-dB prepulse trials to the pulse-alone
      startle amplitude.}
  }
}
\source{
  Data were retrieved and compiled from
  \url{http://dx.doi.org/10.5061/dryad.2rs41}.
}
\references{
  See \link{cfwlab-package} for a list of references.
}
\examples{
# View the first few rows of the phenotype data table.
data(cfw.pheno)
head(cfw.pheno)
}
\keyword{datasets}
