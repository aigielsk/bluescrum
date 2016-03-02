var numeric = require("numeric");

function pinv(A) {
  var z = numeric.svd(A), foo = z.S[0];
  var U = z.U, S = z.S, V = z.V;
  var m = A.length, n = A[0].length, tol = Math.max(m,n)*numeric.epsilon*foo,M = S.length;
  var i,Sinv = new Array(M);
  for(i=M-1;i!==-1;i--) { if(S[i]>tol) Sinv[i] = 1/S[i]; else Sinv[i] = 0; }
  return numeric.dot(numeric.dot(V,numeric.diag(Sinv)),numeric.transpose(U))
}

function findPoints(ms) {
  var A = [];
  var B = [];
  for (i = 0; i < ms.length; ++i) {
    A.push([1, -2 * ms[i][0], -2 * ms[i][1]]);
    B.push(ms[i][2] * ms[i][2] - ms[i][0] * ms[i][0] - ms[i][1] * ms[i][1]);  
  }
  return numeric.dot(pinv(A), B);
}

exports.findPoints = findPoints;
