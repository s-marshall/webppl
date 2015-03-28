////////////////////////////////////////////////////////////////////
// Asynchronous Anytime SMC.
// http://arxiv.org/abs/1407.2864
// bufferSize: queue size
// numParticles: total number of particles to run

'use strict';

var _ = require('underscore');
var util = require('../util.js');
var erp = require('../erp.js');

module.exports = function(env) {

  function copyOneParticle(particle) {
    return {
      continuation: particle.continuation,
      weight: particle.weight,
      completed: particle.completed,
      factorIndex: particle.factorIndex,
      value: particle.value,
      numChildrenToSpawn: 1,
      multiplicity: particle.multiplicity,
      store: _.clone(particle.store)
    };
  }

  function initParticle(s, cont) {
    return {
      continuation: cont,
      weight: 0,
      completed: false,
      factorIndex: undefined,
      value: undefined,
      numChildrenToSpawn: 0,
      multiplicity: 1,
      store: _.clone(s)
    };
  }

  function aSMC(s, k, a, wpplFn, numParticles, bufferSize) {

    this.numParticles = 0;      // K_0 -- initialized here, set in control
    this.bufferSize = bufferSize == undefined ? numParticles : bufferSize; // \rho
    this.initNumParticles = Math.floor(this.bufferSize * (1 / 2));         // \rho_0
    this.exitK = function(s) {return wpplFn(s, env.exit, a);};
    this.store = s;
    this.buffer = [];
    for (var i = 0; i < this.initNumParticles; i++) {
      this.buffer.push(initParticle(this.store, this.exitK));
    }

    this.obsWeights = {};
    this.exitedParticles = [];

    // Move old coroutine out of the way and install this as current handler.
    this.k = k;
    this.oldCoroutine = env.coroutine;
    env.coroutine = this;
    this.oldStore = _.clone(s); // will be reinstated at the end

    // start running
    return this.control(numParticles)
  };

  aSMC.prototype.control = function(numP) {
    console.log("******************* At Control ********************");
    console.log("this.buffer.length =", this.buffer.length);
    // allows for continuation - WIP
    this.numParticles = (numP == undefined) ? this.numParticles : this.numParticles + numP;

    // launch a new particle OR continue an existing one
    var p, launchP;
    var i = Math.floor((this.buffer.length + 1) * Math.random())
    if (i == this.buffer.length) { // generate new particle
      console.log("------------ Generating new particle --------------");
      p = initParticle(this.store, this.exitK);
    } else {                    // launch particle in queue
      launchP = this.buffer[i];
      if (launchP.numChildrenToSpawn > 1) {
        console.log("============== Launching one of many ==============");
        console.log("spawing 1 of " + launchP.numChildrenToSpawn)
        p = copyOneParticle(launchP);
        launchP.numChildrenToSpawn -= 1;
      } else {
        console.log("============== Launching one of one ===============");
        p = launchP;
        this.buffer = util.deleteIndex(this.buffer, i);
      }
    }
    this.activeParticle = p;
    return p.continuation(p.store)
  };

  aSMC.prototype.sample = function(s, cc, a, erp, params) {
    return cc(s, erp.sample(params));
  };

  aSMC.prototype.factor = function(s, cc, a, score) {
    this.activeParticle.weight += score;
    this.activeParticle.continuation = cc;
    this.activeParticle.store = s;
    var fi = this.activeParticle.factorIndex
    var newFI = fi == undefined ? 0 : fi + 1;
    this.activeParticle.factorIndex = newFI;

    var lk = this.obsWeights[newFI];
    if (lk == undefined) {      // 1st particle at observation
      console.log("++++++++++ " + "Particle 1 at factor " + newFI + " ++++++++++");
      var det = {
        wbar: this.activeParticle.weight,
        mnk:  1
      };
      console.log("currWeight = ", this.activeParticle.weight);
      this.obsWeights[newFI] = [det];
      this.activeParticle.numChildrenToSpawn = 1;
    } else {                    // 2nd or greater particle at observation
      console.log("++++++++++ " + "Particle " + (lk.length + 1) + " at factor " + newFI + " ++++++++++");
      var currMultiplicity = this.activeParticle.multiplicity;
      console.log("currMultiplicity = ", currMultiplicity);
      var currWeight = this.activeParticle.weight;
      console.log("currWeight = ", currWeight);
      var denom = lk.length + currMultiplicity; // k - 1 + Ckn
      console.log("lk.length = ", lk.length);
      console.log("denom = ", denom);
      var prevWBar = lk[lk.length-1].wbar;
      console.log("prevWBar = ", prevWBar);
      var wbar = -Math.log(denom) + util.logsumexp([Math.log(lk.length) + prevWBar,
                                                    Math.log(currMultiplicity) + currWeight])
      console.log("wbar = ", wbar);
      if (wbar > 0) throw "Positive weight!!"
      var logRatio = currWeight - wbar;
      console.log("logRatio = ", logRatio);
      var numChildrenAndWeight = [];
      if (logRatio < 0) {
        numChildrenAndWeight = Math.log(Math.random()) < logRatio ?
          [1, wbar] :
          [0, -Infinity];
      } else {
        var totalChildren = 0;
        for (var v = 0; v < lk.length; v++) totalChildren += lk[v].mnk // \sum M^k_n
        console.log("totalChildren = ", totalChildren);
        var minK = Math.min(this.numParticles, lk.length); // min(K_0, k-1)
        var rnk = Math.exp(logRatio);
        var clampedRnk = totalChildren <= minK ? Math.ceil(rnk) : Math.floor(rnk);
        numChildrenAndWeight = [clampedRnk, currWeight - Math.log(clampedRnk)];
      }
      var det2 = {
        wbar: wbar,
        mnk:  numChildrenAndWeight[0]
      };
      console.log("numChildrenAndWeight = ", numChildrenAndWeight);
      this.obsWeights[newFI] = lk.concat([det2])
      if (numChildrenAndWeight[0] > 0) {            // there are children
        if (this.buffer.length < this.bufferSize) { // buffer can be added to
          this.activeParticle.numChildrenToSpawn = numChildrenAndWeight[0];
          this.activeParticle.weight = numChildrenAndWeight[1];
        } else {                  // buffer full, update multiplicty
          this.activeParticle.multiplicity *= numChildrenAndWeight[0];
          this.activeParticle.numChildrenToSpawn = 1;
          this.activeParticle.weight = numChildrenAndWeight[1];
        }
        this.buffer.push(this.activeParticle); // push into buffer
      } else {
        // if there are no children, active particle automatically dropped
        console.log("##### Particle Killed #####")
      }
    }
    return this.control()              // return to control
  };

  aSMC.prototype.exit = function(s, retval) {
    console.log("///// Particle Exited /////")
    this.activeParticle.value = retval;
    this.activeParticle.completed = true;

    // correct weight with multiplicity
    this.activeParticle.weight += Math.log(this.activeParticle.multiplicity)
    this.exitedParticles.push(this.activeParticle)

    if (this.exitedParticles.length < this.numParticles) {
      return this.control()
    } else {
      // Compute marginal distribution from (unweighted) particles
      var hist = {};
      _.each(
        this.exitedParticles,
        function(particle) {
          var k = JSON.stringify(particle.value);
          if (hist[k] === undefined) {
            hist[k] = { prob: 0, val: particle.value };
          }
          hist[k].prob += 1;
        });
      var dist = erp.makeMarginalERP(hist);

      var lastFactorIndex = this.exitedParticles[0].factorIndex;
      console.log("lastFactorIndex = ", lastFactorIndex);
      var lk = this.obsWeights[lastFactorIndex];
      console.log("K_n = ", lk.length);
      console.log("K_0 = ", this.numParticles);
      console.log("Wbar = ", lk[lk.length - 1].wbar);
      dist.normalizationConstant = Math.log(lk.length) // K_n
        - Math.log(this.numParticles)                  // K_0
        + lk[lk.length - 1].wbar;                      // Wbar^k_n

      // allow for continuing pf - WIP
      var ccontrol = this.control.bind(env.coroutine);
      dist.continue = function(numP) {return ccontrol(numP)};

      // Reinstate previous coroutine:
      env.coroutine = this.oldCoroutine;
      // Return from particle filter by calling original continuation:
      return this.k(this.oldStore, dist);
    }
  };

  function asyncPF(s, cc, a, wpplFn, numParticles, bufferSize) {
    return new aSMC(s, cc, a, wpplFn, numParticles, bufferSize);
  }

  return {
    aSMC: asyncPF
  };

};
