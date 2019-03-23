//<!--XXXXXXXXXXXXXXXXXXXXXXXXXX-->
//<!--Copyright 2012 Matt Bognar-->
//<!--XXXXXXXXXXXXXXXXXXXXXXXXXXX-->

//<!--Written by:-->
//<!--Matt Bognar-->
//<!--Department of Statistics and Actuarial Science-->
//<!--University of Iowa-->
//<!--Iowa City, IA 52242-->
//<!--This software may not be re-distributed without the authors consent.-->


function factorial(n) {
	if(n==0 || n==1) {
		return 1;
	}
	else {
		return factorial(n-1)*n;
	}
} 


function logfactorial(n) {
	if(n==0 || n==1) {
		return 0;
		}
	else {
		return logfactorial(n-1) + Math.log(n);
		}
} 


function productRange(a,b) {
  var product=a,i=a;

  while (i++<b) {
    product*=i;
  }

  return product;
}


function combination(n,k) {
  if (n==k || k==0) {
    return 1;
  } 
  else if (n < k) {
  	return 0;
  }
  else {
    k=Math.max(k,n-k);
    return productRange(k+1,n)/productRange(1,n-k);
  }
}


function Fmt(x) { 
	var v

	if(x>=0) { 
		v=''+(x+0.000005) 
		} 
	else { 
		v=''+(x-0.000005) 
		}

	return v.substring(0,v.indexOf('.')+6)
}


function roundNumber(num, dec) 
{
	return(Math.round(num*Math.pow(10,dec)) / Math.pow(10,dec));
}


