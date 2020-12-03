package org.apache.galwat


import org.scalameter.utils.Statistics

import scala.collection.immutable
import scala.collection.immutable.SortedMap
import scala.collection.parallel.CollectionConverters._

/**
 * Node can me considered a male ancestor in Galton-Watson experiment
 *
 * In every new generation it produce 0..n descendant based on Poisson distribution
 *
 * @param lambdaForPoisson The lambda value of the Possion random distrbution generator
 *                         It equals to the average male ancestor for any individual (node)
 */
class Node(val lambdaForPoisson:Double) {


  /**
   * Generates a random number with Poisson distribution
   *
   * Ref : https://stackoverflow.com/questions/1241555/algorithm-to-generate-poisson-and-binomial-random-numbers
   *
   * @param lambda Lambda for Poisson distribution
   * @return
   */
  def nextRandomDescendantsPoisson(lambda:Double) :Int = {
    val l = Math.exp(-lambda)
    var p = 1.0
    var k = 0

    do {
      k += 1
      p *= Math.random
    } while (p > l)
    k-1
  }

  def createChildren() : List[Node] = {
    val numberOfChildren = this.nextRandomDescendantsPoisson(this.lambdaForPoisson)
    List.fill(numberOfChildren)(new Node(this.lambdaForPoisson))
  }
}

/**
 * One Galton-Watson Experiment starting with a seedNode
 *
 * @param maxPopulation The maximum omount of noded that can exist at same time
 *                      This is needed for hardware efficiency too
 *                      One experiment run ends once our lineage got extinct or
 *                      it reaches the max_population
 * @param seedNode The root Node or the Ancestor who lineage is examined
 */
class GaltonWatson(val maxPopulation:Long= 100, val seedNode:Node = new Node(lambdaForPoisson = 1.0)) {
  var livingNodes:immutable.List[Node]= immutable.List[Node]()

  private var _time = 0L
  def time(): Long = _time

  private var _isSeedDominant = false

  /**
   * Has our seedNode became the one and only in the population
   * @return
   */
  def isSeedDominant():Boolean = _isSeedDominant

  def run() :GaltonWatson = {
    // print('.')
    livingNodes = seedNode :: livingNodes
    while(this.livingNodes.nonEmpty && !this.isSeedDominant ) {
      var nextGenNodes = List[Node]()
      for(node <- livingNodes) {
        val children = node.createChildren()
        if (children.size + nextGenNodes.size < maxPopulation)
          nextGenNodes = nextGenNodes ::: children
        else
          this._isSeedDominant = true

      }
      this.livingNodes = nextGenNodes
      _time +=1
    }
    this
  }
}

object GaltonWatson {
  def durationFromMillisToHumanReadable(duration: Long):String = {
    val milliseconds  =  duration % 1000L
    val seconds       = (duration / 1000L) % 60L
    val minutes       = (duration / (1000L*60L)) % 60L
    val hours         = (duration / (1000L*3600L)) % 24L
    val days          = (duration / (1000L*86400L)) % 7L
    val weeks         = (duration / (1000L*604800L)) % 4L
    val months        = (duration / (1000L*2592000L)) % 52L
    val years         = (duration / (1000L*31556952L)) % 10L

    val sb = new scala.collection.mutable.StringBuilder()

    if (years > 0) sb.append(years + " years ")
    if (months > 0) sb.append(months + " months ")
    if (weeks > 0) sb.append(weeks + " weeks ")
    if (days > 0) sb.append(days + " days ")
    if (hours > 0) sb.append(hours + " hours ")
    if (minutes > 0) sb.append(minutes + " min ")
    if (seconds > 0) sb.append(seconds + "s ")
    if (minutes < 1 && hours < 1 && days < 1) {
      if (sb.nonEmpty) sb.append(" ")
      sb.append(milliseconds + "ms")
    }
    sb.toString().trim
  }
  def time[R](block: => R): R = {
    val t0 = System.currentTimeMillis()
    val result = block    // call-by-name
    val t1 = System.currentTimeMillis()
    print("Elapsed time: " + durationFromMillisToHumanReadable(t1-t0))

    result
  }
  def main(args : Array[String]): Unit = {

    time {
      // creating the dimension of runs of our Galton-Watson runs
      // We run n runs in with all of the lambda values described in the cycle
      val dimensions = for(
        lambda <- BigDecimal(1.0) to BigDecimal(1.6) by BigDecimal(0.1);
        n <- 1 to 1000
      ) yield (lambda,n)

      // The required confidence level for the resulting probabilities of the survivals
      val confidence = 0.95

      println(s"Galton-Watson Simulation started with ${dimensions.size} trials")

      val results = dimensions.par
        .map({case (lambda,n)
        => new GaltonWatson(maxPopulation = 1000, seedNode = new Node(lambda.doubleValue)) })
        .map(gw => gw.run())

      // -- calculate and show survival probabilities at various lambdas
      case class ProbabilityWithConfidence(probability:Double,confidence:Double, low:Double, high:Double)
      val survivalProbabilities =  SortedMap[Double,Double]() ++ results.groupBy(gw => gw.seedNode.lambdaForPoisson)
        .map({ case (lambda, gws ) => (lambda,gws.count(gw => gw.isSeedDominant()).toDouble / gws.size)})

      val survivalProbabilitiesConfidence =  SortedMap[Double,ProbabilityWithConfidence]() ++
        results.groupBy(gw => gw.seedNode.lambdaForPoisson)
          .map({ case (lambda, gws ) => (lambda, gws.map(gw => if(gw.isSeedDominant()) 1.0 else 0.0))})
          .map({ case (lambda, zerosAndOnes )
          => (lambda, ProbabilityWithConfidence(probability=zerosAndOnes.reduce(_+_)/zerosAndOnes.size,
            confidence=confidence,
            Statistics.confidenceInterval(List() ++ zerosAndOnes, 1-confidence)._1,
            Statistics.confidenceInterval(List() ++ zerosAndOnes, 1-confidence)._2 ))})



      println("\nSurvival Probabilities Conf Intervals by lambda " )
      survivalProbabilitiesConfidence.map({ case (lambda, ProbabilityWithConfidence(probability, alpha, low, high) )
      => println(s"  - P(survival|lambda=$lambda) = $probability   with ${confidence*100}% conf interval: ["+f"$low%1.2f" + " -> "+ f"$probability%1.2f" +" -> "+f"$high%1.2f" +"]")})

      // -- calculate and show expected extinction time by lambda
      println("\nExpected Extinction time by lambda  ")
      val expExtinctionTime =  SortedMap[Double,Double]()++
        results.groupBy(gw => gw.seedNode.lambdaForPoisson)
          .map({ case (lambda, gws ) => (lambda,gws.filter(! _.isSeedDominant() ).map(gw => gw.time()).reduce(_+_).toDouble / gws.count(! _.isSeedDominant() ))})


      expExtinctionTime.foreach(t=> println(s"  -  E(t(*)|lambda=${t._1}) = ${t._2}"))


    }
  }
}

