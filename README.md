wgs
===

### 说明

capSMART2.0是对`/share/work3/capsmart/pipeline/capSMART/CAPcSMART/capSMART/capSMART_8nt.sjm.nojoin.current.py`

的一个封装，将各程序模块进行独立拆分，并将其放在对应的文件夹中。可用于流程中替代`capSMART_8nt.sjm.nojoin.current.py`.

使用方法和运行环境与`capSMART_8nt.sjm.nojoin.current.py`基本一致。


### 下载和安装:

无需安装，直接源码运行即可：

git clone ssh://git@10.100.14.39:10022/dengyong/capSMART2.0.git

### 基本参数和用法：

python3 /path/wgs/WGS.py -i fastq.txt  -p pair.txt -o  run -c config.txt -L  IDT37.bed

```
 -h, --help            show this help message and exit
  -i LIST, --list LIST  fastq file list
  -b BAMLIST, --bamList BAMLIST
                        bam file list
  -p PAIR, --pair PAIR  pair file list
  -n NUMBERJOB, --numberJob NUMBERJOB
                        the max synchronic job number
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -a, --noanalysis      no analysis
  -c CONFIG, --config CONFIG
                        config file
  -L BED, --bed BED     One genomic intervals over which to operate
  -s STEP, --step STEP  run step snv,cnv,sv,ana
```

详细示例可参考`example`目录中的`run.sh`

生成的`run/.*.job`文件，通过make运行.

### 开发注意事项：

1. 之后大家开发时，脚本中不要出现软件、数据库等的绝对路径，统一从配置文件中读取

    python两种配置文件读取方式：

    方法一、与csmart/QC/同级目录中的脚本可以通过以下命令读取配置文件

    ```
    from utils import *
    from config import config
    conf = Config(config)
    confDict = conf.getdict()
    # 举例获取python路径
    python = confDict['python']
    ```

    方法二、与csmart/QC/scripts同级目录中的脚本可以通过以下命令读取配置文件

    ```
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    script_base_dir = BASE_DIR.split('/csmart/')[0]
    # 将流程根目录添加到PYTHON_PATH中
    sys.path.append(script_base_dir)
    #sys.path.append(os.path.join(script_base_dir,'config')) # python3 需要多添加一个路径到sys.path
    from utils import *
    from config import config

    conf = Config(config)
    confDict = conf.getdict()
    # 举例获取bedtools路径
    bedtools = confDict['bedtools']
    ```

    perl配置文件读取方式：

    方法、在csmart/目录及子目录下的脚本通过以下命令读取配置文件

    ```
    use FindBin '$Bin';
    use Config::IniFiles;
    my $BASE_DIR=(split("/csmart/",$Bin))[0];
    my $cfg = Config::IniFiles->new( -file => "$BASE_DIR/config/config.conf");
    # 举例获取bedtools路径
    my $bedtools ||= $cfg ->val('bin','bedtools');
    ```

    shell.sh配置文件读取方式：

    方法、在csmart/目录及子目录下的脚本通过以下命令读取配置文件
    ```
    BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    BASEDIR=${BASEDIR%%csmart*}
    config=$BASEDIR/config/config.conf
    # __readINI [配置文件路径+名称] [节点名] [键值]
    function __readINI() {
      INIFILE=$1
      SECTION=$2
      ITEM=$3
      _readIni=`awk -F '=' '/\['$SECTION'\]/{a=1}a==1&&$1~/\s?'$ITEM'\s?/&&$1!~/\s?#/{print $2;exit}' $INIFILE`
      echo ${_readIni}
    }
    # 举例获取python路径
    python=( $( __readINI $config bin python ) ) 
    ```

2. 大家需要使用新软件和数据库时，测试时可以安装到自己目录下，测试完成之后，可以联系liuw4318@berryoncology.com拷贝到前边提到的目录下

3. 之后开发的所有脚本，都要统一放到GitLab中，建议该级目录csmart/QC/下只放主脚本，csmart/QC/scripts下放该模块的子脚本


### 更新日志

#### version 4.3.1 
+ 2023-07-24

```
1. 支持52和56 panel
2. 654支持HRD
3. 定制化MRD给MRD_SnvIndel添加DescNew,Abstract,NewMutFunc,并添加DCE3.0注释
```
