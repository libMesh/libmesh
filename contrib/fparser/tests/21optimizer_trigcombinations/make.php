<?php

function CreateTest($f, $file)
{
  $F = $f;
  $C = preg_replace('/(sinh?|cosh?|tanh?|pow|exp)/', 'fp_$1', $f);

  $vars = Array('x');
  #if(strpos($f, '(x)') !== false) $vars[] = 'x';
  #if(strpos($f, '(y)') !== false) $vars[] = 'y';

  print "$F\n";

/*  file_put_contents($file,
    "T=d\n".
    "V=".join(',', $vars)."\n".
    "R=-4,4,0.5\n".
    "F= $F\n".
    "C= $C\n");*/
}

$functions = Array('sin','cos','tan','sinh','cosh','tanh','exp');

for($operator=0; $operator<2; ++$operator)
for($a=0; $a<7; ++$a)
  for($b=$a+1; $b<7; ++$b)
  {
    $f1 = $functions[$a];
    $f2 = $functions[$b];

    for($ae=-2; $ae<=2; ++$ae)
      for($be=-2; $be<=2; ++$be)
      {
        if($be == 0 && $ae == 1) continue; // testing the function alone is not very cool
        if($ae == 0 && $be == 1) continue; // testing the function alone is not very cool

        if($a < 3 && $f2=='exp') continue; // don't bother mixing exp with sin/cos/tan
        if($b < 3 && $f1=='exp') continue; // don't bother mixing exp with sin/cos/tan

        if(!$ae && !$be) continue;

        $afunc = "$f1(x)";
        if($ae == 0) $afunc = "1";
        elseif($ae != 1) $afunc = "pow($afunc,{$ae}.0)";

        $func = $afunc;

        if($be < 0 && $operator==0)
        {
          $func .= " / ";

          $bfunc = "$f2(x)";
          if($be == -1) {}
          else $bfunc = "pow($bfunc,".(-$be).".0)";
          $func .= $bfunc;
        }
        else if($be != 0)
        {
          $func .= ($operator==0 ? " * " : " + ");

          $bfunc = "$f2(x)";
          if($be == 1) {}
          else $bfunc = "pow($bfunc,{$be}.0)";
          $func .= $bfunc;
        }

        static $counter = 0;
        ++$counter;
        $name = $counter;
        #$op = $operator?'add':'mul';
        #$name = sprintf('%s_%sp%d_%sp%d', $op,$f1,$ae,$f2,$be);

        CreateTest($func, $name);
        if(preg_match('/h\(/', $func))
        {
          CreateTest(str_replace('(x)', '(x*x)', $func),
                     ++$counter);
        }
      }
  }
