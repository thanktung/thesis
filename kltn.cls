\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{kltn}

% =====================================================================
% 1. Base Class
% =====================================================================
\LoadClassWithOptions{report}

% =====================================================================
% 2. Ngôn ngữ và mã hóa
% =====================================================================
\RequirePackage[utf8]{vietnam}
\RequirePackage[slovene,english]{babel}
\makeatletter\AtBeginDocument{\let\@elt\relax}\makeatother
% \RequirePackage{csquotes}

% =====================================================================
% 3. Gói cơ bản và toán học
% =====================================================================
\RequirePackage{enumerate}
\RequirePackage{amsmath,amsxtra,amssymb,latexsym,amscd,amsthm}
\RequirePackage{indentfirst,makeidx,setspace}
\RequirePackage[fontsize=13pt]{scrextend}
\setstretch{1.35}

% =====================================================================
% 4. Định dạng đoạn văn
% =====================================================================
\setlength{\parindent}{1cm}
\setlength{\parskip}{6pt}
\setlength{\headsep}{12pt}

% =====================================================================
% 5. Hình ảnh và đồ họa
% =====================================================================
\RequirePackage{graphicx,color,tikz}
\usetikzlibrary{calc}
\RequirePackage{picinpar,floatflt,epic,curves}

% =====================================================================
% 6. Tiêu đề, mục lục và kiểu chữ
% =====================================================================
\RequirePackage{titlesec}
\titleformat{\chapter}{\bfseries \large \center}{CHƯƠNG \thechapter.}{0.3em}{}[]
\titleformat{\section}{\bfseries \large}{ \thesection.}{0.3em}{}[]
\titleformat{\subsection}{\it \bfseries }{ \thesubsection.}{0.3em}{}[]
\titleformat{\subsubsection}{ \it }{ \thesubsubsection.}{0.3em}{}[]
\titlespacing{\chapter}{1em}{0.1em}{2em}
\titlespacing{\section}{-0.5em}{1.5em}{1em}
\titlespacing{\subsection}{-0.5em}{1.3em}{0.5em}
\titlespacing{\subsubsection}{-0.5em}{1.3em}{0.5em}

\RequirePackage[titles]{tocloft}
\setlength\cftchapnumwidth{7em}
\renewcommand\cftchappresnum{CHƯƠNG~}
\renewcommand{\cftchapaftersnum}{.}
\renewcommand{\cfttoctitlefont}{\bfseries \large \hfil}
\renewcommand{\cftloftitlefont}{\bfseries \large \hfil}

% =====================================================================
% 7. Định lý và định nghĩa
% =====================================================================
\RequirePackage{thmtools,tcolorbox,enumitem}
\declaretheoremstyle[
  headfont=\bfseries\scshape,
  notebraces={[}{]},
  bodyfont=\normalfont\itshape,
  headpunct={}
]{myexamplestyle}
\declaretheorem[style=myexamplestyle,name=Định lý,numberwithin=section]{theorem}
\declaretheorem[style=myexamplestyle]{definition}
\newtheorem{defi}[theorem]{Định nghĩa}
\newtheorem{lemma}[theorem]{Bổ đề}

% =====================================================================
% 8. Thuật toán
% =====================================================================
\RequirePackage[chapter]{algorithm}
\RequirePackage{algorithmic}
\floatname{algorithm}{\hfil Thuật toán}

% =====================================================================
% 9. Header/Footer
% =====================================================================
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyheadoffset{0cm}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
\fancyfoot[C]{\thepage}
\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyfoot[C]{\thepage}%
}

% =====================================================================
% 10. Hyperref và Bookmarks
% =====================================================================
\RequirePackage[
  hidelinks,
  unicode,
  bookmarksopen,
  bookmarksnumbered
]{hyperref}
\RequirePackage{bookmark,eso-pic,calc}

% =====================================================================
% 11. Công thức và hiển thị
% =====================================================================
\allowdisplaybreaks
\makeatletter
\g@addto@macro\normalsize{%
  \setlength\abovedisplayskip{6pt}
  \setlength\belowdisplayskip{6pt}
  \setlength\abovedisplayshortskip{6pt}
  \setlength\belowdisplayshortskip{6pt}%
}
\makeatother

% =====================================================================
% 12. Hình ảnh và bảng biểu
% =====================================================================
\RequirePackage{booktabs}
\RequirePackage[format=plain,labelfont={bf,it},textfont=it]{caption}

% =====================================================================
% 13. Phông chữ đặc biệt và lệnh custom
% =====================================================================
\newcommand\HRule{\rule{\textwidth}{1pt}}
\renewcommand\qedsymbol{\ensuremath{\square}}
% \font\tua=ugqbo8v at 30pt % đã thay thế bằng font Unicode chuẩn
\newcommand*{\MyDef}{\mathrm{def}}
\newcommand*{\eqdefU}{\ensuremath{\mathop{\overset{\MyDef}{=}}}}
\newcommand*{\myeq}{\mathop{\overset{\MyDef}{\resizebox{\widthof{\eqdefU}}{\heightof{=}}{=}}}}

% =====================================================================
% 14. Trang, mục lục, phụ lục
% =====================================================================
\RequirePackage[a4paper,left=3.0cm, right=2cm,top=2.5cm,bottom=3cm]{geometry}
\RequirePackage{pdfpages}

\newcommand\VKnumRoman{\pagenumbering{roman}\setcounter{page}{1}}
\newcommand\VKbatDaudanhSo{\newpage\pagenumbering{arabic}\setcounter{page}{1}}
\newcommand\VKngatTrang{\newpage\renewcommand{\cleardoublepage}{}\renewcommand{\clearpage}{}}
\newcommand\VKdanhSoPhuLuc{\newpage\renewcommand*{\thepage}{P\arabic{page}}\setcounter{page}{1}}

\newcommand\VKmucLuc{\newpage\renewcommand*\contentsname{MỤC LỤC}\tableofcontents}
\newcommand\VKdanhMucHinhVe{\renewcommand*{\listfigurename}{\bfseries DANH MỤC CÁC HÌNH VẼ}{\let\oldnumberline\numberline\renewcommand{\numberline}{Hình~\oldnumberline}\listoffigures}}
\newcommand\VKdanhMucBangBieu{\renewcommand*{\listtablename}{\bfseries DANH MỤC CÁC BẢNG BIỂU}{\let\oldnumberline\numberline\renewcommand{\numberline}{Bảng~\oldnumberline}\listoftables}}

\makeatletter
\g@addto@macro\appendix{\addtocontents{toc}{\protect\renewcommand{\protect\cftchappresnum}{PHỤ LỤC~}}}
\makeatother

% =====================================================================
% 15. Tài liệu tham khảo
% =====================================================================
\RequirePackage[backend=biber, style=numeric-comp, sorting=nyvt, maxbibnames=120, block=none, autolang=other]{biblatex}

\DefineBibliographyStrings{slovene}{
  andothers = {và cộng sự},
  pages     = {tr.},
  and       = {và},
  editor    = {Biên tập},
  editors   = {Các biên tập},
  in        = {trong},
  chapter   = {Chương},
  edition   = {lần xuất bản},
  urlseen   = {truy cập}
}

\newcommand\VKtaiLieuThamKhao{\chapter*{TÀI LIỆU THAM KHẢO}\addcontentsline{toc}{chapter}{{\bf TÀI LIỆU THAM KHẢO}}\titleformat{\chapter}{\bfseries  \large}{CHƯƠNG \thechapter.}{}{}[]\titlespacing{\chapter}{0.1em}{0.1em}{1.1em}
\printbibliography[notkeyword={Vietnam},heading=subbibliography,title={Tiếng Anh:}]
\printbibliography[keyword={Vietnam},heading=subbibliography,title={Tiếng Việt:}]}

% =====================================================================
% 16. Danh mục khác
% =====================================================================
\newcommand\VKdanhMucDinhLy{\renewcommand{\listtheoremname}{DANH MỤC CÁC ĐỊNH LÝ}\addcontentsline{toc}{section}{{\bf DANH MỤC CÁC ĐỊNH LÝ\rm }}\listoftheorems [ignoreall,show={theorem}]}
\newcommand\VKdanhMucDinhNghia{\renewcommand{\listtheoremname}{DANH MỤC CÁC ĐỊNH NGHĨA}\addcontentsline{toc}{section}{{\bf DANH MỤC CÁC ĐỊNH NGHĨA\rm }}\listoftheorems [ignoreall,show={defi}]}
\newcommand\VKdanhMucListings{\addcontentsline{toc}{section}{{\bf DANH MỤC VIẾT TẮT\rm }}}

% =====================================================================
% 17. Gói mở rộng cho bảng, hình, mã nguồn
% =====================================================================
\RequirePackage{commath,multirow,adjustbox,longtable,float,tabularx,ragged2e}

% =====================================================================
% 18. Môi trường listings cho mã nguồn
% =====================================================================
\RequirePackage{listings}
\RequirePackage{color}
\definecolor{mygreen}{RGB}{28,172,0}
\definecolor{mylilas}{RGB}{170,55,241}

\lstset{
  language=Matlab,
  basicstyle=\footnotesize\ttfamily,
  breaklines=true,
  morekeywords={matlab2tikz},
  keywordstyle=\color{blue},
  identifierstyle=\color{black},
  stringstyle=\color{mylilas},
  commentstyle=\color{mygreen},
  showstringspaces=false,
  numbers=left,
  numberstyle={\tiny \color{black}},
  numbersep=9pt,
  emph=[1]{for,end,break},
  emphstyle=[1]\color{red}
}

% =====================================================================
% 19. Kiểu cột bảng p{width} can trái, giữa, phải
% =====================================================================
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}
\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}
\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}

% =====================================================================
% 20. Unicode ký tự đặc biệt
% =====================================================================
\DeclareUnicodeCharacter{2212}{-}
\DeclareUnicodeCharacter{0301}{+}
\DeclareUnicodeCharacter{0308}{+}
\DeclareUnicodeCharacter{00B1}{+}
