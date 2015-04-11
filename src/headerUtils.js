'use strict';

var util = require('./util.js');

module.exports = function(env) {

  function display(s, k, a, x) {
    return k(s, console.log(x));
  }

  // Caching for a wppl function f.
  //
  // Caution!: if `f` isn't deterministic weird stuff can happen.
  // caching is across all uses of `f`, even in different execution paths.
  function cache(s, k, a, f) {
    var c = util.initHashMap();
    var cf = function(s, k, a) {
      var args = Array.prototype.slice.call(arguments, 3);
      var lookup = c.get(args);
      if (lookup) {
        return k(s, lookup);
      } else {
        var newk = function(s, r) {
          c.set(args, r);
          return k(s, r);
        };
        return f.apply(this, [s, newk, a].concat(args));
      }
    };
    return k(s, cf);
  }

  function apply(s, k, a, wpplFn, args) {
    return wpplFn.apply(global, [s, k, a].concat(args));
  }

  return {
    display: display,
    cache: cache,
    apply: apply
  };

};
